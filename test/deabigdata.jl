# Tests for Big Data Radial DEA Models
@testset "BigData RadialDEAModel" begin

    ## Test Radial DEA Models with FLS Book data
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    # Input oriented CRS
    deaio = deabigdata(X, Y, orient = :Input, rts = :CRS)

    @test typeof(deaio) == RadialDEAModel

    @test nobs(deaio) == 11
    @test ninputs(deaio) == 2
    @test noutputs(deaio) == 1
    @test efficiency(deaio) ≈ [1.0000000000;
                               0.6222896791;
                               0.8198562444;
                               1.0000000000;
                               0.3103709311;
                               0.5555555556;
                               1.0000000000;
                               0.7576690896;
                               0.8201058201;
                               0.4905660377;
                               1.0000000000]
    @test convert(Matrix, peers(deaio)) ≈
    [1.000000000  0  0 0.0000000000  0  0 0.00000000000  0  0   0   0;
     0.000000000  0  0 0.4249783174  0  0 0.10928013877  0  0   0   0;
     1.134321653  0  0 0.4380053908  0  0 0.00000000000  0  0   0   0;
     0.000000000  0  0 1.0000000000  0  0 0.00000000000  0  0   0   0;
     0.000000000  0  0 0.2573807721  0  0 0.04844814534  0  0   0   0;
     0.000000000  0  0 0.0000000000  0  0 0.33333333333  0  0   0   0;
     0.000000000  0  0 0.0000000000  0  0 1.00000000000  0  0   0   0;
     0.000000000  0  0 1.0348650979  0  0 0.11457435013  0  0   0   0;
     0.000000000  0  0 0.0000000000  0  0 1.14814814815  0  0   0   0;
     0.000000000  0  0 0.4905660377  0  0 0.49056603774  0  0   0   0;
     0.000000000  0  0 0.0000000000  0  0 0.00000000000  0  0   0   1.000000000] atol = 1e-8
     @test slacks(deaio, :X) ≈ [0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                4.444444444   0;
                                0.000000000   0;
                                0.000000000   0;
                                1.640211640   0;
                                0.000000000   0;
                                0.000000000   4]
    @test slacks(deaio, :Y) ≈ zeros(11)

    @test efficiency(deabigdata(targets(deaio, :X), targets(deaio, :Y), orient = :Input, rts = :CRS, slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaio, :X), targets(deaio, :Y))) ≈ zeros(11) atol=1e-14

    # Otuput oriented CRS
    deaoo = deabigdata(X, Y, orient = :Output, rts = :CRS)

    @test nobs(deaoo) == 11
    @test ninputs(deaoo) == 2
    @test noutputs(deaoo) == 1
    @test efficiency(deaoo) ≈ [1.0000000000;
                               1.606968641;
                               1.219726027;
                               1.0000000000;
                               3.221951220;
                               1.800000000;
                               1.0000000000;
                               1.319837398;
                               1.219354839;
                               2.038461538;
                               1.0000000000]
    @test convert(Matrix, peers(deaoo)) ≈
    [1.000000000  0  0 0.0000000000  0  0 0.00000000000  0  0   0   0;
     0.000000000  0  0 0.6829268293  0  0 0.1756097561   0  0   0   0;
     1.383561644  0  0 0.5342465753  0  0 0.00000000000  0  0   0   0;
     0.000000000  0  0 1.0000000000  0  0 0.00000000000  0  0   0   0;
     0.000000000  0  0 0.8292682927  0  0 0.1560975610   0  0   0   0;
     0.000000000  0  0 0.0000000000  0  0 0.6000000000   0  0   0   0;
     0.000000000  0  0 0.0000000000  0  0 1.00000000000  0  0   0   0;
     0.000000000  0  0 1.3658536585  0  0 0.1512195122   0  0   0   0;
     0.000000000  0  0 0.0000000000  0  0 1.4000000000   0  0   0   0;
     0.000000000  0  0 1.0000000000  0  0 1.0000000000   0  0   0   0;
     1.000000000  0  0 0.0000000000  0  0 0.00000000000  0  0   0   0]
     @test slacks(deaoo, :X) ≈ [0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                8   0;
                                0.000000000   0;
                                0.000000000   0;
                                2   0;
                                0.000000000   0;
                                0.000000000   4]
    @test slacks(deaoo, :Y) ≈ zeros(11)

    @test efficiency(deabigdata(targets(deaoo, :X), targets(deaoo, :Y), orient = :Output, rts = :CRS, slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaoo, :X), targets(deaoo, :Y))) ≈ zeros(11) atol=1e-10

    # Input oriented VRS
    deaiovrs = deabigdata(X, Y, orient = :Input, rts = :VRS)

    @test nobs(deaiovrs) == 11
    @test ninputs(deaiovrs) == 2
    @test noutputs(deaiovrs) == 1
    @test efficiency(deaiovrs) ≈ [1.0000000000;
                                  0.8699861687;
                                  1.0000000000;
                                  1.0000000000;
                                  0.7116402116;
                                  1.0000000000;
                                  1.0000000000;
                                  1.0000000000;
                                  1.0000000000;
                                  0.4931209269;
                                  1.0000000000]
    @test convert(Matrix, peers(deaiovrs)) ≈
    [1.000000000    0  0 0.0000000000  0  0.00000000000 0.00000000000  0  0   0   0;
     0.52558782849  0  0 0.0000000000  0  0.2842323651  0.1901798064   0  0   0   0;
     0.000000000    0  1 0.0000000000  0  0.00000000000 0.00000000000  0  0   0   0;
     0.000000000    0  0 1.0000000000  0  0.00000000000 0.00000000000  0  0   0   0;
     0.56613756614  0  0 0.0000000000  0  0.4338624339  0.00000000000  0  0   0   0;
     0.000000000    0  0 0.0000000000  0  1.00000000000 0.00000000000  0  0   0   0;
     0.000000000    0  0 0.0000000000  0  0.00000000000 1.00000000000  0  0   0   0;
     0.000000000    0  0 0.0000000000  0  0.00000000000 0.00000000000  1  0   0   0;
     0.000000000    0  0 0.0000000000  0  0.00000000000 0.00000000000  0  1   0   0;
     0.03711078928  0  0 0.4433381608  0  0.00000000000 0.5195510500   0  0   0   0;
     0.000000000    0  0 0.0000000000  0  0.00000000000 0.00000000000  0  0   0   1.000000000]
     @test slacks(deaiovrs, :X) ≈ [0.000000000   0;
                                   0.000000000   0;
                                   0.000000000   0;
                                   0.000000000   0;
                                   0.000000000   0;
                                   0   0;
                                   0.000000000   0;
                                   0.000000000   0;
                                   0   0;
                                   0.000000000   0;
                                   0.000000000   4]
    @test slacks(deaiovrs, :Y) ≈ [0.000000000;
                                  0.000000000;
                                  0.000000000;
                                  0.000000000;
                                  2.698412698;
                                  0.000000000;
                                  0.000000000;
                                  0.000000000;
                                  0.000000000;
                                  0.000000000;
                                  0.000000000]

    @test efficiency(deabigdata(targets(deaiovrs, :X), targets(deaiovrs, :Y), orient = :Input, rts = :VRS, slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaiovrs, :X), targets(deaiovrs, :Y))) ≈ zeros(11) atol=1e-12

    # Output oriented VRS
    deaoovrs = deabigdata(X, Y, orient = :Output, rts = :VRS)

    @test nobs(deaoovrs) == 11
    @test ninputs(deaoovrs) == 2
    @test noutputs(deaoovrs) == 1
    @test efficiency(deaoovrs) ≈ [1.0000000000;
                                  1.507518797;
                                  1.0000000000;
                                  1.0000000000;
                                  3.203947368;
                                  1.000000000;
                                  1.0000000000;
                                  1.000000000;
                                  1.000000000;
                                  1.192307692;
                                  1.0000000000]
    @test convert(Matrix, peers(deaoovrs)) ≈
    [1.000000000   0  0 0.0000000000  0  0 0.00000000000  0  0   0   0;
     0.38157894737 0  0 0.1710526316  0  0 0.4473684211   0  0   0   0;
     0.000000000   0  1 0.0000000000  0  0 0.00000000000  0  0   0   0;
     0.000000000   0  0 1.0000000000  0  0 0.00000000000  0  0   0   0;
     0.03947368421 0  0 0.7763157895  0  0 0.1842105263   0  0   0   0;
     0.000000000   0  0 0.0000000000  0  1 0.00000000000  0  0   0   0;
     0.000000000   0  0 0.0000000000  0  0 1.00000000000  0  0   0   0;
     0.000000000   0  0 0.0000000000  0  0 0.00000000000  1  0   0   0;
     0.000000000   0  0 0.0000000000  0  0 0.00000000000  0  1   0   0;
     0.000000000   0  0 0.0000000000  0  0 0.00000000000  0  1   0   0;
     1.000000000   0  0 0.0000000000  0  0 0.00000000000  0  0   0   0]
    @test slacks(deaoovrs, :X) ≈ [0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                5   11;
                                0.000000000   4]
    @test slacks(deaoovrs, :Y) ≈ zeros(11) atol=1e-10

    @test efficiency(deabigdata(targets(deaoovrs, :X), targets(deaoovrs, :Y), orient = :Output, rts = :VRS, slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaoovrs, :X), targets(deaoovrs, :Y))) ≈ zeros(11) atol=1e-12

    # Test no slacks
    deaionoslack = deabigdata(X, Y, slack = false)
    @test efficiency(deaionoslack) == efficiency(deaio)
    @test isempty(slacks(deaionoslack, :X)) == 1
    @test isempty(slacks(deaionoslack, :Y)) == 1

    @test efficiency(deabigdata(targets(deaionoslack, :X), targets(deaionoslack, :Y), slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaionoslack, :X), targets(deaionoslack, :Y))) != zeros(11) # Different as there is no slacks in first model

    # Print
    show(IOBuffer(), deaio)
    show(IOBuffer(), deaionoslack)

    # Test errors
    @test_throws DimensionMismatch deabigdata([1; 2 ; 3], [4 ; 5]) #  Different number of observations
    @test_throws ArgumentError deabigdata([1; 2; 3], [4; 5; 6], orient = :Error) # Invalid orientation
    @test_throws ArgumentError deabigdata([1; 2; 3], [4; 5; 6], rts = :Error) # Invalid returns to scale

    # ------------------
    # Test if no exteriors
    # ------------------
    Xnoext = [1 1; 1.5 1; 2 1]
    Ynoext = [2 2; 1.5 1.5; 1 0.5]

    deanoext = dea(Xnoext, Ynoext, orient = :Input)
    deabignoext = deabigdata(Xnoext, Ynoext, orient = :Input)

    @test efficiency(deanoext) ≈ efficiency(deabignoext)
    @test slacks(deanoext, :X) ≈ slacks(deabignoext, :X)
    @test slacks(deanoext, :X) ≈ slacks(deabignoext, :X)
    @test targets(deanoext, :X) ≈ targets(deabignoext, :X)
    @test targets(deanoext, :Y) ≈ targets(deabignoext, :Y)
    @test peersmatrix(deanoext) ≈ peersmatrix(deabignoext)

    # ------------------
    # Test with random data
    # ------------------
    rng = StableRNG(1234567)
    X = rand(Uniform(10, 20), 500, 6)
    Y = rand(Uniform(10, 20), 500, 4)

    rdea = dea(X, Y, progress = false)
    rdeabig = deabigdata(X, Y, progress = false)

    @test efficiency(rdeabig) ≈ efficiency(rdea)
    @test slacks(rdeabig, :X) ≈ slacks(rdea, :X)
    @test slacks(rdeabig, :Y) ≈ slacks(rdea, :Y)
    @test targets(rdeabig, :X) ≈ targets(rdea, :X)
    @test targets(rdeabig, :Y) ≈ targets(rdea, :Y)
    @test peersmatrix(rdeabig) ≈ peersmatrix(rdea)

    # ------------------
    # Test Vector and Matrix inputs and outputs
    # ------------------
    # Tests against results in R

    # Inputs is Matrix, Outputs is Vector
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deabigdata(X, Y, orient = :Input)) ≈ [1; 1; 1; 0.6; 0.4; 1; 0.6666666667; 0.625]

    # Inputs is Vector, Output is Matrix
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]

    @test efficiency(deabigdata(X, Y, orient = :Output)) ≈ [1; 1; 1; 1.555555556; 2.333333333; 1; 1.272727273; 1.6]

    # Inputs is Vector, Output is Vector
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]

    @test efficiency(deabigdata(X, Y, orient = :Input)) ≈ [0.4; 1; 0.8; 0.6; 0.4; 0.4; 0.5142857143; 0.2]

end
