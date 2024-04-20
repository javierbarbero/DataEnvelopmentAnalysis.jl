# Tests for Directional DEA Models
@testset "DirectionalDEAModel" begin

    ## Test Directional DF DEA Models with FLS Book data
    # Test against results in R
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    # Observed CRS
    deaddfobs = deaddf(X, Y, Gx = X, Gy = Y, rts = :CRS)

    @test typeof(deaddfobs) == DirectionalDEAModel

    @test nobs(deaddfobs) == 11
    @test ninputs(deaddfobs) == 2
    @test noutputs(deaddfobs) == 1
    @test efficiency(deaddfobs) ≈ [0.00000000000;
                               0.23282544774;
                               0.09898790422;
                               0.00000000000;
                               0.52628538417;
                               0.28571428571;
                               0.00000000000;
                               0.13787061050;
                               0.09883720930;
                               0.34177215190;
                               0.00000000000]
    @test slacks(deaddfobs, :X) ≈ [0.000000000   0;
                               0.000000000   0;
                               0.000000000   0;
                               0.000000000   0;
                               0.000000000   0;
                               5.714285714   0;
                               0.000000000   0;
                               0.000000000   0;
                               1.802325581   0;
                               0.000000000   0;
                               0.000000000   4]
    @test slacks(deaddfobs, :Y) ≈ zeros(11)

    @test efficiency(deaddf(X, Y, Gx = :Observed, Gy = :Observed, rts = :CRS)) == efficiency(deaddfobs)

    @test efficiency(deaddf(targets(deaddfobs, :X), targets(deaddfobs, :Y), Gx = :Observed, Gy = :Observed, rts = :CRS, slack = false)) ≈ zeros(11) atol=1e-15
    @test efficiency(deaadd(targets(deaddfobs, :X), targets(deaddfobs, :Y))) ≈ zeros(11) atol=1e-13

    # Observed VRS
    deaddfobsvrs = deaddf(X, Y, Gx = X, Gy = Y, rts = :VRS)

    @test nobs(deaddfobsvrs) == 11
    @test ninputs(deaddfobsvrs) == 2
    @test noutputs(deaddfobsvrs) == 1
    @test efficiency(deaddfobsvrs) ≈ [0.00000000000;
                              0.1076130509;
                              0.0000000000;
                              0.0000000000;
                              0.2883597884;
                              0.0000000000;
                              0.00000000000;
                              0.0000000000;
                              0.0000000000;
                              0.1821192053;
                              0.00000000000]
    @test slacks(deaddfobsvrs, :X) ≈ [0.000000000   0;
                              0.000000000   0;
                              0.000000000   0;
                              0.000000000   0;
                              0.000000000   0;
                              0   0;
                              0.000000000   0;
                              0.000000000   0;
                              0   0;
                              0.000000000   4.32781457;
                              0.000000000   4]
    @test slacks(deaddfobsvrs, :Y) ≈ [0.000000000;
                             0.000000000;
                             0.000000000;
                             0.000000000;
                             0.3915343915;
                             0.000000000;
                             0.000000000;
                             0.000000000;
                             0.000000000;
                             0.000000000;
                             0.000000000]

    @test efficiency(deaddf(X, Y, Gx = :Observed, Gy = :Observed, rts = :VRS)) == efficiency(deaddfobsvrs)

    @test efficiency(deaddf(targets(deaddfobsvrs, :X), targets(deaddfobsvrs, :Y), Gx = :Observed, Gy = :Observed, rts = :VRS, slack = false)) ≈ zeros(11) atol=1e-15
    @test efficiency(deaadd(targets(deaddfobsvrs, :X), targets(deaddfobsvrs, :Y))) ≈ zeros(11) atol=1e-12

    # Ones CRS
    deaddfones = deaddf(X, Y, Gx = ones(size(X)), Gy = ones(size(Y)), rts = :CRS)

    @test nobs(deaddfones) == 11
    @test ninputs(deaddfones) == 2
    @test noutputs(deaddfones) == 1
    @test efficiency(deaddfones) ≈ [0.00000000000;
                             3.219963031;
                             2.121693122;
                             0.00000000000;
                             6.735674677;
                             1.945945946;
                             0.00000000000;
                             3.635859519;
                             1.837837838;
                             10.231053604;
                             0.00000000000]
    @test slacks(deaddfones, :X) ≈ [0.000000000   0;
                             0.000000000   0;
                             0.000000000   0;
                             0.000000000   0;
                             0.000000000   0;
                             10.918918919   0;
                             0.000000000   0;
                             0.000000000   0;
                             4.756756757    0;
                             0.000000000   0;
                             0.000000000   4]
    @test slacks(deaddfones, :Y) ≈ zeros(11)

    @test efficiency(deaddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :CRS)) == efficiency(deaddfones)

    # Ones VRS
    deaddfonesvrs = deaddf(X, Y, Gx = ones(size(X)), Gy = ones(size(Y)), rts = :VRS)

    @test nobs(deaddfonesvrs) == 11
    @test ninputs(deaddfonesvrs) == 2
    @test noutputs(deaddfonesvrs) == 1
    @test efficiency(deaddfonesvrs) ≈ [0.000000000;
                            1.418867925;
                            0.000000000;
                            0.000000000;
                            4.067924528;
                            0.000000000;
                            0.000000000;
                            0.000000000;
                            0.000000000;
                            5.000000000;
                            0.000000000]
    @test slacks(deaddfonesvrs, :X) ≈ [0.000000000   0;
                            0.000000000   0;
                            0.000000000   0;
                            0.000000000   0;
                            0.000000000   0;
                            0   0;
                            0.000000000   0;
                            0.000000000   0;
                            0   0;
                            0.000000000   6;
                            0.000000000   4]
    @test slacks(deaddfonesvrs, :Y) ≈ zeros(11) atol = 1e-14

    @test efficiency(deaddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :VRS)) == efficiency(deaddfonesvrs)

    # Only X CRS
    deaddfonlyX = deaddf(X, Y, Gx =  X, Gy = zeros(size(Y)), rts = :CRS)

    @test nobs(deaddfonlyX) == 11
    @test ninputs(deaddfonlyX) == 2
    @test noutputs(deaddfonlyX) == 1
    @test efficiency(deaddfonlyX) ≈ [0.0000000000;
                            0.3777103209;
                            0.1801437556;
                            0.0000000000;
                            0.6896290689;
                            0.4444444444;
                            0.0000000000;
                            0.2423309104;
                            0.1798941799;
                            0.5094339623;
                            0.0000000000]
    @test slacks(deaddfonlyX, :X) ≈ [0.000000000   0;
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
    @test slacks(deaddfonlyX, :Y) ≈ zeros(11)

    @test efficiency(deaddf(X, Y, Gx = :Observed, Gy = :Zeros, rts = :CRS)) == efficiency(deaddfonlyX)

    # Only X VRS
    deaddfonlyXvrs = deaddf(X, Y, Gx = X, Gy = zeros(size(Y)), rts = :VRS)

    @test nobs(deaddfonlyXvrs) == 11
    @test ninputs(deaddfonlyXvrs) == 2
    @test noutputs(deaddfonlyXvrs) == 1
    @test efficiency(deaddfonlyXvrs) ≈ [0.0000000000;
                            0.1300138313;
                            0.0000000000;
                            0.0000000000;
                            0.2883597884;
                            0.0000000000;
                            0.0000000000;
                            0.0000000000;
                            0.0000000000;
                            0.5068790731;
                            0.0000000000]
    @test slacks(deaddfonlyXvrs, :X) ≈ [0.000000000   0;
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
    @test slacks(deaddfonlyXvrs, :Y) ≈ [0.000000000;
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

    @test efficiency(deaddf(X, Y, Gx = :Observed, Gy = :Zeros, rts = :VRS)) == efficiency(deaddfonlyXvrs)

    # Only Y CRS
    deaddfonlyY = deaddf(X, Y, Gx = zeros(size(X)), Gy = Y, rts = :CRS)

    @test nobs(deaddfonlyY) == 11
    @test ninputs(deaddfonlyY) == 2
    @test noutputs(deaddfonlyY) == 1
    @test efficiency(deaddfonlyY) ≈ [0.0000000000;
                            0.6069686411;
                            0.2197260274;
                            0.0000000000;
                            2.2219512195;
                            0.8000000000;
                            0.0000000000;
                            0.3198373984;
                            0.2193548387;
                            1.0384615385;
                            0.0000000000]
    @test slacks(deaddfonlyY, :X) ≈ [0.000000000   0;
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
    @test slacks(deaddfonlyY, :Y) ≈ zeros(11)

    @test efficiency(deaddf(X, Y, Gx = :Zeros, Gy = :Observed, rts = :CRS)) == efficiency(deaddfonlyY)

    # Only Y VRS
    deaddfonlyYvrs = deaddf(X, Y, Gx = zeros(size(X)), Gy = Y, rts = :VRS)

    @test nobs(deaddfonlyYvrs) == 11
    @test ninputs(deaddfonlyYvrs) == 2
    @test noutputs(deaddfonlyYvrs) == 1
    @test efficiency(deaddfonlyYvrs) ≈ [0.0000000000;
                            0.5075187970;
                            0.0000000000;
                            0.0000000000;
                            2.2039473684;
                            0.0000000000;
                            0.0000000000;
                            0.0000000000;
                            0.0000000000;
                            0.1923076923;
                            0.0000000000]
    @test slacks(deaddfonlyYvrs, :X) ≈ [0.000000000   0;
                            0.000000000   0;
                            0.000000000   0;
                            0.000000000   0;
                            0.000000000   0;
                            0   0;
                            0.000000000   0;
                            0.000000000   0;
                            0   0;
                            5   11;
                            0.000000000   4]
    @test slacks(deaddfonlyYvrs, :Y) ≈ zeros(11) atol = 1e-14

    @test efficiency(deaddf(X, Y, Gx = :Zeros, Gy = :Observed, rts = :VRS)) == efficiency(deaddfonlyYvrs)

    # Mean CRS
    deaddfmean = deaddf(X, Y, Gx = :Mean, Gy = :Mean, rts = :CRS)

    @test efficiency(deaddfmean) ≈ [0; 0.1713509018; 0.1082533683; 0; 0.3584401184; 0.1148158887; 0; 0.1934829069; 0.1084372282; 0.5444473258; 0]

    # Mean VRS
    deaddfmeanvrs = deaddf(X, Y, Gx = :Mean, Gy = :Mean, rts = :VRS)

    @test efficiency(deaddfmeanvrs) ≈ [0; 0.08056567388; 0; 0;0.23098350118; 0; 0; 0; 0; 0.24886877828; 0]

    # Observed FDH
    deaddffdh = deaddf(X, Y, Gx = X, Gy = Y, rts = :FDH)

    @test efficiency(deaddffdh) ≈ [0.0000000000;
                                0.0000000000;
                                0.0000000000;
                                0.0000000000;
                                0.1111111111;
                                0.0000000000;
                                0.0000000000;
                                0.0000000000;
                                0.0000000000;
                                0.1200000000;
                                0.0000000000]
    @test slacks(deaddffdh, :X) ≈ [0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0.4444444444444464;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                0.000000000   0;
                                9.96   0;
                                0.000000000   4]
    @test slacks(deaddffdh, :Y) ≈ [0.0; 0.0; 0.0; 0.0; 5.11111111111112; 0.0; 0.0; 0.0; 0.0; 0.879999999999999; 0.0]

    @test efficiency(deaddf(targets(deaddffdh, :X), targets(deaddffdh, :Y), Gx = :Observed, Gy = :Observed, rts = :FDH, slack = false)) ≈ zeros(11) atol=1e-15
    @test efficiency(deaadd(targets(deaddffdh, :X), targets(deaddffdh, :Y), rts = :FDH)) ≈ zeros(11) atol=1e-13

    # Test no slacks
    deaddfnoslack = deaddf(X, Y, Gx = X, Gy = Y, slack = false)
    @test efficiency(deaddfnoslack) == efficiency(deaddfobs)
    @test isempty(slacks(deaddfnoslack, :X)) == 1
    @test isempty(slacks(deaddfnoslack, :Y)) == 1

    @test efficiency(dea(targets(deaddfnoslack, :X), targets(deaddfnoslack, :Y), slack = false)) ≈ ones(11)
    @test efficiency(deaadd(targets(deaddfnoslack, :X), targets(deaddfnoslack, :Y))) != zeros(11) # Different as there is no slacks in first model

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deaddfobs_ref_eff = zeros(size(X, 1))

    deaddfobsvs_ref_eff = zeros(size(X, 1))

    deaddfvrs_ref_slackX = zeros(size(X))
    deaddfvrs_ref_slackY = zeros(size(Y))

    Gx = X[:,:]
    Gy = Y[:,:]

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]
        Gxeval = Gx[i:i,:]
        Gxeval = Gxeval[:,:]
        Gyeval = Gy[i:i,:]
        Gyeval = Gyeval[:,:]

        deaddfobs_ref_eff[i] = efficiency(deaddf(Xeval, Yeval, Gx = Gxeval, Gy = Gyeval, rts = :CRS, Xref = Xref, Yref = Yref))[1]

        deaddfobsvs_ref_eff[i] = efficiency(deaddf(Xeval, Yeval, Gx = Gxeval, Gy = Gyeval, rts = :VRS, Xref = Xref, Yref = Yref))[1]

        deaddfvrs_ref_slackX[i,:] = slacks(deaddf(Xeval, Yeval, Gx = Gxeval, Gy = Gyeval, rts = :VRS, Xref = Xref, Yref = Yref), :X)
        deaddfvrs_ref_slackY[i,:] = slacks(deaddf(Xeval, Yeval, Gx = Gxeval, Gy = Gyeval, rts = :VRS, Xref = Xref, Yref = Yref), :Y)
    end

    @test deaddfobs_ref_eff ≈ efficiency(deaddfobs)

    @test deaddfobsvs_ref_eff ≈ efficiency(deaddfobsvrs)

    @test deaddfvrs_ref_slackX ≈ slacks(deaddfobsvrs, :X) atol=1e-14
    @test deaddfvrs_ref_slackY ≈ slacks(deaddfobsvrs, :Y) atol=1e-15

    # Print
    show(IOBuffer(), deaddfobs)
    show(IOBuffer(), deaddfnoslack)

    # Test errors
    @test_throws DimensionMismatch deaddf([1; 2 ; 3], [4 ; 5], Gx = [1; 2 ; 3], Gy = [4 ; 5]) #  Different number of observations
    @test_throws DimensionMismatch deaddf([1; 2], [4 ; 5], Gx = [1; 2], Gy = [4 ; 5], Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws DimensionMismatch deaddf([1 1; 2 2], [4 4; 5 5], Gx = [1 1; 2 2], Gy = [4 4; 5 5], Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws DimensionMismatch deaddf([1 1; 2 2], [4 4; 5 5], Gx = [1 1; 2 2], Gy = [4 4; 5 5], Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ArgumentError deaddf([1; 2; 3], [4; 5; 6], Gx = [1; 2; 3], Gy = [4; 5; 6], rts = :Error) # Invalid returns to scale
    @test_throws DimensionMismatch deaddf([1 1; 2 2; 3 3], [4; 5; 6], Gx = [1 1 1; 2 2 2; 3 3 3], Gy = [4; 5; 6]) # Different size of inputs direction
    @test_throws DimensionMismatch deaddf([1; 2; 3], [4 4; 5 5; 6 6], Gx = [1; 2; 3], Gy = [4 4 4; 5 5 5; 6 6 6]) # Different size of inputs direction
    @test_throws ArgumentError deaddf([1; 2; 3], [1; 2; 3], Gx = :Error, Gy = :Ones) # Invalid inuts direction
    @test_throws ArgumentError deaddf([1; 2; 3], [1; 2; 3], Gx = :Ones, Gy = :Error) # Invalid outputs direction

    # ------------------
    # Test Vector and Matrix inputs and outputs
    # ------------------
    # Tests against results in R

    # Inputs is Matrix, Outputs is Vector
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deaddf(X, Y, Gx = X, Gy = Y)) ≈ [0; 0; 0; 0.25; 0.4285714286; 0; 0.2; 0.2307692308]

    # Inputs is Vector, Output is Matrix
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]

    @test efficiency(deaddf(X, Y, Gx = X, Gy = Y)) ≈ [0; 0; 0; 0.2173913043; 0.4; 0; 0.12; 0.2307692308]

    # Inputs is Vector, Output is Vector
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]

    @test efficiency(deaddf(X, Y, Gx = X, Gy = Y)) ≈ [0.4285714286; 0; 0.1111111111; 0.25; 0.4285714286; 0.4285714286; 0.3207547170; 0.6666666667]

    # ------------------
    # Test model with zero Directions
    # ------------------   
    deamddfzeros = deaddf(X, Y, Gx = zeros(size(X)), Gy = zeros(size(Y)), rts = :VRS)
    @test efficiency(deamddfzeros) == zeros(8)


end
