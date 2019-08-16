# Tests for Radial DEA Models
@testset "DirectionalDEAModel" begin

    ## Test Directional DF DEA Models with FLS Book data
    # Test against results with R
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    # Observed CRS
    deaddfobs = deaddf(X, Y, X, Y, rts = :CRS)

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

    # Observed VRS
    deaddfobsvrs = deaddf(X, Y, X, Y, rts = :VRS)

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

    # Ones CRS
    deaddfones = deaddf(X, Y, ones(size(X)), ones(size(Y)), rts = :CRS)

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

    # Ones VRS
    deaddfonesvrs = deaddf(X, Y, ones(size(X)), ones(size(Y)), rts = :VRS)

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

    # Only X CRS
    deaddfonlyX = deaddf(X, Y, X, zeros(size(Y)), rts = :CRS)

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

    # Only X VRS
    deaddfonlyXvrs = deaddf(X, Y, X, zeros(size(Y)), rts = :VRS)

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

    # Only Y CRS
    deaddfonlyY = deaddf(X, Y, zeros(size(X)), Y, rts = :CRS)

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

    # Only Y VRS
    deaddfonlyYvrs = deaddf(X, Y, zeros(size(X)), Y, rts = :VRS)

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

     ## Test if one-by-one DEA using evaluation and reference sets match initial results
     deaddfobs_ref_eff = zeros(size(X, 1))

     deaddfobsvs_ref_eff = zeros(size(X, 1))

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

         deaddfobs_ref_eff[i] = efficiency(deaddf(Xeval, Yeval, Gxeval, Gyeval, rts = :CRS, Xref = Xref, Yref = Yref))[1]

         deaddfobsvs_ref_eff[i] = efficiency(deaddf(Xeval, Yeval, Gxeval, Gyeval, rts = :VRS, Xref = Xref, Yref = Yref))[1]
     end

     @test deaddfobs_ref_eff ≈ efficiency(deaddfobs)

     @test deaddfobsvs_ref_eff ≈ efficiency(deaddfobsvrs)

end
