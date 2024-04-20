# Tests for Radial DEA Models
@testset "RadialDEAModel" begin

    ## Test Radial DEA Models with FLS Book data
    # Test against results in R
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    # Input oriented CRS
    deaio = dea(X, Y, orient = :Input, rts = :CRS)

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

    @test efficiency(dea(targets(deaio, :X), targets(deaio, :Y), orient = :Input, rts = :CRS, slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaio, :X), targets(deaio, :Y))) ≈ zeros(11) atol=1e-14

    @test peersmatrix(deaio) == deaio.lambda

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

    @test efficiency(dea(targets(deaoo, :X), targets(deaoo, :Y), orient = :Output, rts = :CRS, slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaoo, :X), targets(deaoo, :Y))) ≈ zeros(11) atol=1e-10

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

    @test efficiency(dea(targets(deaiovrs, :X), targets(deaiovrs, :Y), orient = :Input, rts = :VRS, slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaiovrs, :X), targets(deaiovrs, :Y))) ≈ zeros(11) atol=1e-12

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
    @test slacks(deaoovrs, :Y) ≈ zeros(11) atol=1e-10

    @test efficiency(dea(targets(deaoovrs, :X), targets(deaoovrs, :Y), orient = :Output, rts = :VRS, slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaoovrs, :X), targets(deaoovrs, :Y))) ≈ zeros(11) atol=1e-12

    # Input oriented FDH
    deafdh = dea(X, Y, orient = :Input, rts = :FDH)

    @test typeof(deafdh) == RadialDEAModel

    @test nobs(deafdh) == 11
    @test ninputs(deafdh) == 2
    @test noutputs(deafdh) == 1
    @test efficiency(deafdh) ≈ [1.0000000000;
                                1.0000000000;
                                1.0000000000;
                                1.0000000000;
                                0.8888888889;
                                1.0000000000;
                                1.0000000000;
                                1.0000000000;
                                1.0000000000;
                                0.5952380952;
                                1.0000000000]
    @test convert(Matrix, peers(deafdh)) ≈
        [1  0  0  0  0  0  0  0  0   0   0;
        0  1  0  0  0  0  0  0  0   0   0;
        0  0  1  0  0  0  0  0  0   0   0;
        0  0  0  1  0  0  0  0  0   0   0;
        0  1  0  0  0  0  0  0  0   0   0;
        0  0  0  0  0  1  0  0  0   0   0;
        0  0  0  0  0  0  1  0  0   0   0;
        0  0  0  0  0  0  0  1  0   0   0;
        0  0  0  0  0  0  0  0  1   0   0;
        0  0  0  0  0  0  1  0  0   0   0;
        0  0  0  0  0  0  0  0  0   0   1]
     @test slacks(deafdh, :X) ≈ [0.0 0.0;
                                 0.0 0.0;
                                 0.0 0.0;
                                 0.0 0.0;
                                 0.0 0.44444444444444287;
                                 0.0 0.0;
                                 0.0 0.0;
                                 0.0 0.0;
                                 0.0 0.0;
                                 3.552713678800501e-15 4.880952380952383;
                                 0.0 4.0]
    @test slacks(deafdh, :Y) ≈ [0.0; 0.0; 0.0; 0.0; 6.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.0];

    @test efficiency(dea(targets(deafdh, :X), targets(deafdh, :Y), orient = :Input, rts = :FDH, slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deafdh, :X), targets(deafdh, :Y), rts = :FDH)) ≈ zeros(11) atol=1e-14

    @test peersmatrix(deafdh) == deafdh.lambda

    ## Test no slacks
    deaionoslack = dea(X, Y, slack = false)
    @test efficiency(deaionoslack) == efficiency(deaio)
    @test isempty(slacks(deaionoslack, :X)) == 1
    @test isempty(slacks(deaionoslack, :Y)) == 1

    @test efficiency(dea(targets(deaionoslack, :X), targets(deaionoslack, :Y), slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaionoslack, :X), targets(deaionoslack, :Y))) != zeros(11) # Different as there is no slacks in first model

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
    @test_throws DimensionMismatch dea([1; 2 ; 3], [4 ; 5]) #  Different number of observations
    @test_throws DimensionMismatch dea([1; 2], [4 ; 5], Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws DimensionMismatch dea([1 1; 2 2], [4 4; 5 5], Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws DimensionMismatch dea([1 1; 2 2], [4 4; 5 5], Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ArgumentError dea([1; 2; 3], [4; 5; 6], orient = :Error) # Invalid orientation
    @test_throws ArgumentError dea([1; 2; 3], [4; 5; 6], rts = :Error) # Invalid returns to scale

    @test_throws ArgumentError dea([1; 2 ; 3], [4 ; 5; 6], disposX = :Error)  # Invalid inputs disposability
    @test_throws ArgumentError dea([1; 2 ; 3], [4 ; 5; 6], disposY = :Error)  # Invalid outputs disposability

    @test_throws ArgumentError targets(deaio, :Error) # Invalid target
    @test_throws ArgumentError slacks(deaio, :Error) # Invalid slacks
    @test_throws ArgumentError rts(deaio) # rts only for DEA models in multiplier form

    # ------------------
    # Weak Disposability Tests
    # ------------------

    X = [1; 2; 3; 2; 4]
    Y = [2; 3; 4; 1; 3]

    deaioStrong = dea(X, Y, orient = :Input, rts = :VRS)
    @test efficiency(deaioStrong ) ≈ [1.0; 1.0; 1.0; 0.5; 0.5]
    @test slacks(deaioStrong, :X) ≈ [0; 0; 0; 0; 0] atol=1e-15
    @test slacks(deaioStrong, :Y) ≈ [0; 0; 0; 1; 0] atol=1e-15

    deaioWeak = dea(X, Y, orient = :Input, rts = :VRS, disposY = :Weak)
    @test efficiency(deaioWeak ) ≈ [1.0; 1.0; 1.0; 1.0; 0.5]
    @test slacks(deaioWeak, :X) ≈ [0; 0; 0; 0; 0] atol=1e-15
    @test slacks(deaioWeak, :Y) ≈ [0; 0; 0; 0; 0] atol=1e-15

    deaooStrong = dea(X, Y, orient = :Output, rts = :VRS)
    @test efficiency(deaooStrong ) ≈ [1.0; 1.0; 1.0; 3.0; 1.3333333333333333]
    @test slacks(deaooStrong, :X) ≈ [0; 0; 0; 0; 1] atol=1e-15
    @test slacks(deaooStrong, :Y) ≈ [0; 0; 0; 0; 0] atol=1e-15

    deaooWeak = dea(X, Y, orient = :Output, rts = :VRS, disposX = :Weak)
    @test efficiency(deaooWeak ) ≈ [1.0; 1.0; 1.0; 3.0; 1.0]
    @test slacks(deaooWeak, :X) ≈ [0; 0; 0; 0; 0] atol=1e-14
    @test slacks(deaooWeak, :Y) ≈ [0; 0; 0; 0; 0] atol=1e-14

    # Test if weak disposability in inputs in the Input oriented model works
    # In this example, same result as stron disposability
    @test efficiency(dea(X, Y, orient = :Input, disposX = :Strong)) ==
        efficiency(dea(X, Y, orient = :Input, disposX = :Weak))

    # Test if weak disposability in outputs in the Output oriented model works
    # In this example, same result as strong disposability
    @test efficiency(dea(X, Y, orient = :Output, disposY = :Strong)) ==
        efficiency(dea(X, Y, orient = :Output, disposY = :Weak))

    # ------------------
    # DMU names
    # ------------------

    X = [1; 2; 3; 2; 4]
    Y = [2; 3; 4; 1; 3]

    @test names(dea(X, Y)) == ["1"; "2"; "3"; "4"; "5"]
    @test names(dea(X, Y, names = ["A", "B", "C", "D", "E"])) == ["A", "B", "C", "D", "E"]

    logs, value = Test.collect_test_logs() do
        names(dea(X, Y, names = ["A", "B", "C", "D"]))
    end
    @test occursin("Length of names lower than number of observations", string(logs))
    @test value == ["A", "B", "C", "D", "5"]

    logs, value = Test.collect_test_logs() do
        names(dea(X, Y, names = ["A", "B", "C", "D", "E", "F"]))
    end
    @test occursin("Length of names greater than number of observations", string(logs))
    @test value == ["A", "B", "C", "D", "E"]

    # ------------------
    # Test Vector and Matrix inputs and outputs
    # ------------------
    # Tests against results in R

    # Inputs is Matrix, Outputs is Vector
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(dea(X, Y, orient = :Input)) ≈ [1; 1; 1; 0.6; 0.4; 1; 0.6666666667; 0.625]

    # Inputs is Vector, Output is Matrix
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]

    @test efficiency(dea(X, Y, orient = :Output)) ≈ [1; 1; 1; 1.555555556; 2.333333333; 1; 1.272727273; 1.6]

    # Inputs is Vector, Output is Vector
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]

    @test efficiency(dea(X, Y, orient = :Input)) ≈ [0.4; 1; 0.8; 0.6; 0.4; 0.4; 0.5142857143; 0.2]

    # ------------------
    # Run with different progress meter options
    # ------------------
    # Just run the code for code coverage as the progress meter does not affect the result
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    dea(X, Y, progress = true)
    dea(X, Y, progress = false)

end
