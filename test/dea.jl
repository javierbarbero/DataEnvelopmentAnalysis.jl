# Tests for Radial DEA Models
@testset "RadialDEAModel" begin

    ## Test Radial DEA Models with FLS Book data
    # Test against results with R
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    # Input oriented CRS
    deaio = dea(X, Y, orient = :Input, rts = :CRS)

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
     0.000000000  0  0 0.0000000000  0  0 0.00000000000  0  0   0   1.000000000]
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

    # Otuput oriented CRS
    deaoo = dea(X, Y, orient = :Output, rts = :CRS)

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

    # Input oriented VRS
    deaiovrs = dea(X, Y, orient = :Input, rts = :VRS)

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

    # Output oriented VRS
    deaoovrs = dea(X, Y, orient = :Output, rts = :VRS)

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
    @test slacks(deaoovrs, :Y) ≈ zeros(11) atol=1e-15
    
    # Test no slacks
    deaionoslack = dea(X, Y, slack = false)
    @test efficiency(deaionoslack) == efficiency(deaio)
    @test isempty(slacks(deaionoslack, :X)) == 1
    @test isempty(slacks(deaionoslack, :Y)) == 1

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deaio_ref_eff = zeros(size(X, 1))
    deaoo_ref_eff = zeros(size(X, 1))

    deaiovrs_ref_eff = zeros(size(X, 1))
    deaoovrs_ref_eff = zeros(size(X, 1))

    deaiovrs_ref_slackX = zeros(size(X))
    deaiovrs_ref_slackY = zeros(size(Y))

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]

        deaio_ref_eff[i] = efficiency(dea(Xeval, Yeval, orient = :Input, rts = :CRS, Xref = Xref, Yref = Yref))[1]
        deaoo_ref_eff[i] = efficiency(dea(Xeval, Yeval, orient = :Output, rts = :CRS, Xref = Xref, Yref = Yref))[1]

        deaiovrs_ref_eff[i] = efficiency(dea(Xeval, Yeval, orient = :Input, rts = :VRS, Xref = Xref, Yref = Yref))[1]
        deaoovrs_ref_eff[i] = efficiency(dea(Xeval, Yeval, orient = :Output, rts = :VRS, Xref = Xref, Yref = Yref))[1]

        deaiovrs_ref_slackX[i,:] = slacks(dea(Xeval, Yeval, orient = :Input, rts = :VRS, Xref = Xref, Yref = Yref), :X)
        deaiovrs_ref_slackY[i,:] = slacks(dea(Xeval, Yeval, orient = :Input, rts = :VRS, Xref = Xref, Yref = Yref), :Y)
    end

    @test deaio_ref_eff ≈ efficiency(deaio)
    @test deaoo_ref_eff ≈ efficiency(deaoo)

    @test deaiovrs_ref_eff ≈ efficiency(deaiovrs)
    @test deaoovrs_ref_eff ≈ efficiency(deaoovrs)

    @test deaiovrs_ref_slackX ≈ slacks(deaiovrs, :X) atol=1e-14
    @test deaiovrs_ref_slackY ≈ slacks(deaiovrs, :Y) atol=1e-15

    # Print
    show(IOBuffer(), deaio)
    show(IOBuffer(), deaionoslack)

    # Test errors
    @test_throws ErrorException dea([1; 2 ; 3], [4 ; 5]) #  Different number of observations
    @test_throws ErrorException dea([1; 2], [4 ; 5], Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws ErrorException dea([1 1; 2 2], [4 4; 5 5], Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws ErrorException dea([1 1; 2 2], [4 4; 5 5], Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ErrorException dea([1; 2; 3], [4; 5; 6], orient = :Error) # Invalid orientation
    @test_throws ErrorException dea([1; 2; 3], [4; 5; 6], rts = :Error) # Invalid returns to scale

end
