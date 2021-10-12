# Tests for Directional Multiplier DEA Models
@testset "DirectionalMultiplierDEAModel" begin

    ## Test Directional DF DEA Models with FLS Book data
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    # Observed CRS
    deaddfmobs = deaddfm(X, Y, Gx = X, Gy = Y, rts = :CRS)

    @test typeof(deaddfmobs) == DirectionalMultiplierDEAModel

    @test nobs(deaddfmobs) == 11
    @test ninputs(deaddfmobs) == 2
    @test noutputs(deaddfmobs) == 1
    @test efficiency(deaddfmobs) ≈ [
        0.00000000000;
        0.23282544774;
        0.09898790422;
        0.00000000000;
        0.52628538417;
        0.28571428571;
        0.00000000000;
        0.13787061050;
        0.09883720930;
        0.34177215190;
        0.00000000000] atol = 1e-5
    @test multipliers(deaddfmobs, :X) ≈ [
        0.0450913    0.0211187
        0.0193798    0.0255279
        0.0195014    0.00913355
        0.0136023    0.0179174
        0.0209417    0.0275852
        0.0          0.107143
        0.0130985    0.0172538
        0.0101633    0.0133875
        0.0          0.0392442
        0.00895338   0.0117938
        0.1         -1.11022e-16
    ] atol = 1e-5
    @test multipliers(deaddfmobs, :Y) ≈ [
        0.04166666666666667
        0.027399091152098366
        0.0180202419155764
        0.01923076923076923
        0.02960716348931253
        0.0396825396825397
        0.018518518518518517
        0.014368823158337424
        0.014534883720930238
        0.012658227848101266
        0.04166666666666667
    ] atol = 1e-5

    @test efficiency(deaddfm(X, Y, Gx = :Observed, Gy = :Observed, rts = :CRS)) == efficiency(deaddfmobs)

    @test efficiency(deaddfm(targets(deaddfmobs, :X), targets(deaddfmobs, :Y), Gx = :Observed, Gy = :Observed, rts = :CRS)) ≈ zeros(11) atol = 1e-5
    @test rts(deaddfmobs) ≈ zeros(11)

    # Observed VRS
    deaddfmobsvrs = deaddfm(X, Y, Gx = X, Gy = Y, rts = :VRS)

    @test nobs(deaddfmobsvrs) == 11
    @test ninputs(deaddfmobsvrs) == 2
    @test noutputs(deaddfmobsvrs) == 1
    @test efficiency(deaddfmobsvrs) ≈ [
        0.00000000000;
        0.1076130509;
        0.0000000000;
        0.0000000000;
        0.2883597884;
        0.0000000000;
        0.00000000000;
        0.0000000000;
        0.0000000000;
        0.1821192053;
        0.00000000000] atol = 1e-5
    @test multipliers(deaddfmobsvrs, :X) ≈ [
        0.0450913    0.0211187
        0.0163137    0.0472238
        0.0264496    0.000339098
        0.0136023    0.0179174
        0.0185185    0.047619
        0.0          0.125
        0.0130985    0.0172538
        0.00627716   0.00482859
        0.0          0.0222222
        0.00331126   0.0
        0.1         -1.11022e-16
    ] atol = 1e-5
    @test multipliers(deaddfmobsvrs, :Y) ≈ [
        0.04166666666666667
        0.012306811677160837
        0.02271956595456087
        0.01923076923076923
        0.0
        0.027777777777777773
        0.018518518518518517
        0.024142926122646062
        0.022222222222222206
        0.033112582781456984
        0.04166666666666667
    ] atol = 1e-5

    @test efficiency(deaddfm(X, Y, Gx = :Observed, Gy = :Observed, rts = :VRS)) == efficiency(deaddfmobsvrs)

    @test efficiency(deaddfm(targets(deaddfmobsvrs, :X), targets(deaddfmobsvrs, :Y), Gx = :Observed, Gy = :Observed, rts = :VRS)) ≈ zeros(11) atol = 1e-5
    @test rts(deaddfmobsvrs) ≈ [
        0.0
        -0.5477962220950202
         0.1359782977280435
         0.0
        -0.7116402116402119
        -0.5000000000000003
         0.0
         0.448575567358764
         0.37777777777777694
         0.9039735099337763
         0.0
    ] atol = 1e-5

    # Ones CRS
    deaddfmones = deaddfm(X, Y, Gx = ones(size(X)), Gy = ones(size(Y)), rts = :CRS)

    @test nobs(deaddfmones) == 11
    @test ninputs(deaddfmones) == 2
    @test noutputs(deaddfmones) == 1
    @test efficiency(deaddfmones) ≈ [
        0.00000000000;
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
    @test multipliers(deaddfmones, :X) ≈ [
        0.705882  0.0
        0.268022  0.35305
        0.417989  0.195767
        0.417989  0.195767
        0.268022  0.35305
        0.0       0.72973
        0.268022  0.35305
        0.268022  0.35305
        0.0       0.72973
        0.268022  0.35305
        0.705882  0.0
    ] atol = 1e-5
    @test multipliers(deaddfmones, :Y) ≈ [
        0.29411764705882354
        0.3789279112754159
        0.3862433862433863
        0.3862433862433863
        0.3789279112754159
        0.27027027027027034
        0.3789279112754159
        0.3789279112754159
        0.27027027027027034
        0.3789279112754159
        0.29411764705882354
    ] atol = 1e-5

    @test efficiency(deaddfm(X, Y, Gx = :Ones, Gy = :Ones, rts = :CRS)) == efficiency(deaddfmones)
    @test rts(deaddfmobs) ≈ zeros(11)

    # Ones VRS
    deaddfmonesvrs = deaddfm(X, Y, Gx = ones(size(X)), Gy = ones(size(Y)), rts = :VRS)

    @test nobs(deaddfmonesvrs) == 11
    @test ninputs(deaddfmonesvrs) == 2
    @test noutputs(deaddfmonesvrs) == 1
    @test efficiency(deaddfmonesvrs) ≈ [
        0.000000000;
        1.418867925;
        0.000000000;
        0.000000000;
        4.067924528;
        0.000000000;
        0.000000000;
        0.000000000;
        0.000000000;
        5.000000000;
        0.000000000] atol = 1e-5
    @test multipliers(deaddfmonesvrs, :X) ≈ [
        0.705882   0.0
        0.215094   0.622642
        0.534247   0.00684932
        0.417989   0.195767
        0.215094   0.622642
        0.0        0.818182
        0.268022   0.35305
        0.178082   0.136986
        0.0        0.5
        0.0909091  0.0
        0.705882   0.0
    ] atol = 1e-5
    @test multipliers(deaddfmonesvrs, :Y) ≈ [
        0.29411764705882354
        0.16226415094339627
        0.4589041095890411
        0.3862433862433863
        0.16226415094339627
        0.18181818181818193
        0.3789279112754159
        0.6849315068493145
        0.4999999999999999
        0.9090909090909107
        0.29411764705882354
    ] atol = 1e-5

    @test efficiency(deaddfm(X, Y, Gx = :Ones, Gy = :Ones, rts = :VRS)) == efficiency(deaddfmonesvrs)
    @test rts(deaddfmonesvrs) ≈ [
        0.0
        -7.222641509433959
         2.746575342465754
         0.0
        -7.222641509433959
        -3.2727272727272703
         0.0
        12.72602739726025
         8.499999999999993
        24.81818181818191
         0.0
    ] atol = 1e-5

    # Only X CRS
    deaddfmonlyX = deaddfm(X, Y, Gx =  X, Gy = zeros(size(Y)), rts = :CRS)

    @test nobs(deaddfmonlyX) == 11
    @test ninputs(deaddfmonlyX) == 2
    @test noutputs(deaddfmonlyX) == 1
    @test efficiency(deaddfmonlyX) ≈ [
        0.0000000000;
        0.3777103209;
        0.1801437556;
        0.0000000000;
        0.6896290689;
        0.4444444444;
        0.0000000000;
        0.2423309104;
        0.1798941799;
        0.5094339623;
        0.0000000000] atol = 1e-5
    @test multipliers(deaddfmonlyX, :X) ≈ [
        0.0901826  0.0422374
        0.0314397  0.0414137
        0.0354897  0.0166217
        0.0416228  0.0194942
        0.0274413  0.0361469
        0.0        0.166667
        0.0261969  0.0345077
        0.0178637  0.0235309
        0.0        0.0714286
        0.0133456  0.0175794
        0.2        2.61229e-17
    ] atol = 1e-5
    @test multipliers(deaddfmonlyX, :Y) ≈ [
        0.08333333333333336
        0.04444926279271466
        0.032794249775381853
        0.038461538461538464
        0.03879636638909917
        0.06172839506172836
        0.03703703703703704
        0.02525563631883701
        0.026455026455026464
        0.018867924528301886
        0.0833333333333333
    ] atol = 1e-5

    @test efficiency(deaddfm(X, Y, Gx = :Observed, Gy = :Zeros, rts = :CRS)) == efficiency(deaddfmonlyX)

    # Only X VRS
    deaddfmonlyXvrs = deaddfm(X, Y, Gx = X, Gy = zeros(size(Y)), rts = :VRS)

    @test nobs(deaddfmonlyXvrs) == 11
    @test ninputs(deaddfmonlyXvrs) == 2
    @test noutputs(deaddfmonlyXvrs) == 1
    @test efficiency(deaddfmonlyXvrs) ≈ [
        0.0000000000;
        0.1300138313;
        0.0000000000;
        0.0000000000;
        0.2883597884;
        0.0000000000;
        0.0000000000;
        0.0000000000;
        0.0000000000;
        0.5068790731;
        0.0000000000] atol = 1e-5
    @test multipliers(deaddfmonlyXvrs, :X) ≈ [
        0.0901826  0.0422374
        0.0197095  0.0570539
        0.0612245  0.000784929
        0.0416228  0.0194942
        0.0185185  0.047619
        0.0        0.166667
        0.0261969  0.0345077
        0.0227671  0.0175131
        0.0        0.0714286
        0.013034   0.0181028
        0.2        2.61229e-17
    ] atol = 1e-5
    @test multipliers(deaddfmonlyXvrs, :Y) ≈ [
        0.08333333333333336
        0.014868603042876911
        0.05259026687598115
        0.038461538461538464
        0.0
        0.03703703703703698
        0.03703703703703704
        0.08756567425569173
        0.07142857142857138
        0.013758146270818254
        0.0833333333333333
    ] atol = 1e-5

    @test efficiency(deaddfm(X, Y, Gx = :Observed, Gy = :Zeros, rts = :VRS)) == efficiency(deaddfmonlyXvrs)

    # Only Y CRS
    deaddfmonlyY = deaddfm(X, Y, Gx = zeros(size(X)), Gy = Y, rts = :CRS)

    @test nobs(deaddfmonlyY) == 11
    @test ninputs(deaddfmonlyY) == 2
    @test noutputs(deaddfmonlyY) == 1
    @test efficiency(deaddfmonlyY) ≈ [
        0.0000000000;
        0.6069686411;
        0.2197260274;
        0.0000000000;
        2.2219512195;
        0.8000000000;
        0.0000000000;
        0.3198373984;
        0.2193548387;
        1.0384615385;
        0.0000000000] atol = 1e-5
    @test multipliers(deaddfmonlyY, :X) ≈ [
        0.0901826  0.0422374
        0.0505226  0.0665505
        0.0432877  0.020274
        0.0272045  0.0358349
        0.0884146  0.116463
        0.0        0.3
        0.0261969  0.0345077
        0.0235772  0.0310569
        0.0        0.0870968
        0.0272045  0.0358349
        0.2        0.0
    ] atol = 1e-5
    @test multipliers(deaddfmonlyY, :Y) ≈ [
        0.08333333333333329
        0.07142857142857142
        0.04000000000000001
        0.03846153846153846
        0.12499999999999999
        0.11111111111111112
        0.037037037037037035
        0.03333333333333334
        0.03225806451612903
        0.03846153846153846
        0.08333333333333318
    ] atol = 1e-5

    @test efficiency(deaddfm(X, Y, Gx = :Zeros, Gy = :Observed, rts = :CRS)) == efficiency(deaddfmonlyY)

    # Only Y VRS
    deaddfmonlyYvrs = deaddfm(X, Y, Gx = zeros(size(X)), Gy = Y, rts = :VRS)

    @test nobs(deaddfmonlyYvrs) == 11
    @test ninputs(deaddfmonlyYvrs) == 2
    @test noutputs(deaddfmonlyYvrs) == 1
    @test efficiency(deaddfmonlyYvrs) ≈ [
        0.0000000000;
        0.5075187970;
        0.0000000000;
        0.0000000000;
        2.2039473684;
        0.0000000000;
        0.0000000000;
        0.0000000000;
        0.0000000000;
        0.1923076923;
        0.0000000000] atol = 1e-5
    @test multipliers(deaddfmonlyYvrs, :X) ≈ [
        0.0901826   0.0422374
        0.0676692   0.093985
        0.0465672   0.000597015
        0.0272045   0.0358349
        0.118421    0.164474
        0.0         0.5
        0.0261969   0.0345077
        0.00866667  0.00666667
        0.0         0.0322581
        0.0         0.0
        0.2         0.0
    ] atol = 1e-5
    @test multipliers(deaddfmonlyYvrs, :Y) ≈ [
        0.08333333333333329
        0.07142857142857144
        0.040000000000000015
        0.03846153846153846
        0.125
        0.11111111111111109
        0.037037037037037035
        0.03333333333333333
        0.03225806451612903
        0.038461538461538464
        0.08333333333333318
    ] atol = 1e-5

    @test efficiency(deaddfm(X, Y, Gx = :Zeros, Gy = :Observed, rts = :VRS)) == efficiency(deaddfmonlyYvrs)

    # Mean CRS
    deaddfmmean = deaddfm(X, Y, Gx = :Mean, Gy = :Mean, rts = :CRS)

    @test efficiency(deaddfmmean) ≈ [0; 0.1713509018; 0.1082533683; 0; 0.3584401184; 0.1148158887; 0; 0.1934829069; 0.1084372282; 0.5444473258; 0] atol = 1e-5

    # Mean VRS
    deaddfmmeanvrs = deaddfm(X, Y, Gx = :Mean, Gy = :Mean, rts = :VRS)

    @test efficiency(deaddfmmeanvrs) ≈ [0; 0.08056567388; 0; 0;0.23098350118; 0; 0; 0; 0; 0.24886877828; 0] atol = 1e-5

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deamddf_ref_eff = zeros(size(X, 1))

    deamddfvrs_ref_eff = zeros(size(X, 1))

    deamddfvrs_ref_multX = zeros(size(X))
    deamddfvrs_ref_multY = zeros(size(Y))

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

        deamddf_ref_eff[i] = efficiency(deaddfm(Xeval, Yeval, Gx = Gxeval, Gy = Gyeval, rts = :CRS, Xref = Xref, Yref = Yref))[1]

        deamddfvrs_ref_eff[i] = efficiency(deaddfm(Xeval, Yeval, Gx = Gxeval, Gy = Gyeval, rts = :VRS, Xref = Xref, Yref = Yref))[1]

        deamddfvrs_ref_multX[i,:] = multipliers(deaddfm(Xeval, Yeval, Gx = Gxeval, Gy = Gyeval, rts = :VRS, Xref = Xref, Yref = Yref), :X)
        deamddfvrs_ref_multY[i,:] = multipliers(deaddfm(Xeval, Yeval, Gx = Gxeval, Gy = Gyeval, rts = :VRS, Xref = Xref, Yref = Yref), :Y)
    end

    @test deamddf_ref_eff ≈ efficiency(deaddfmobs)

    @test deamddfvrs_ref_eff ≈ efficiency(deaddfmobsvrs)

    @test deamddfvrs_ref_multX ≈ multipliers(deaddfmobsvrs, :X) atol = 1e-5
    @test deamddfvrs_ref_multY ≈ multipliers(deaddfmobsvrs, :Y) atol = 1e-5

    # Print
    show(IOBuffer(), deaddfmobs)

    # Test errors
    @test_throws DimensionMismatch deaddfm([1; 2 ; 3], [4 ; 5], Gx = [1; 2 ; 3], Gy = [4 ; 5]) #  Different number of observations
    @test_throws DimensionMismatch deaddfm([1; 2], [4 ; 5], Gx = [1; 2], Gy = [4 ; 5], Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws DimensionMismatch deaddfm([1 1; 2 2], [4 4; 5 5], Gx = [1 1; 2 2], Gy = [4 4; 5 5], Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws DimensionMismatch deaddfm([1 1; 2 2], [4 4; 5 5], Gx = [1 1; 2 2], Gy = [4 4; 5 5], Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ArgumentError deaddfm([1; 2; 3], [4; 5; 6], Gx = [1; 2; 3], Gy = [4; 5; 6], rts = :Error) # Invalid returns to scale
    @test_throws DimensionMismatch deaddfm([1 1; 2 2; 3 3], [4; 5; 6], Gx = [1 1 1; 2 2 2; 3 3 3], Gy = [4; 5; 6]) # Different size of inputs direction
    @test_throws DimensionMismatch deaddfm([1; 2; 3], [4 4; 5 5; 6 6], Gx = [1; 2; 3], Gy = [4 4 4; 5 5 5; 6 6 6]) # Different size of inputs direction
    @test_throws ArgumentError deaddfm([1; 2; 3], [1; 2; 3], Gx = :Error, Gy = :Ones) # Invalid inuts direction
    @test_throws ArgumentError deaddfm([1; 2; 3], [1; 2; 3], Gx = :Ones, Gy = :Error) # Invalid outputs direction

    # ------------------
    # Test Vector and Matrix inputs and outputs
    # ------------------
    # Tests against results in R

    # Inputs is Matrix, Outputs is Vector
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deaddfm(X, Y, Gx = X, Gy = Y)) ≈ [0; 0; 0; 0.25; 0.4285714286; 0; 0.2; 0.2307692308] atol = 1e-5

    # Inputs is Vector, Output is Matrix
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]

    @test efficiency(deaddfm(X, Y, Gx = X, Gy = Y)) ≈ [0; 0; 0; 0.2173913043; 0.4; 0; 0.12; 0.2307692308] atol = 1e-5

    # Inputs is Vector, Output is Vector
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]

    @test efficiency(deaddfm(X, Y, Gx = X, Gy = Y)) ≈ [0.4285714286; 0; 0.1111111111; 0.25; 0.4285714286; 0.4285714286; 0.3207547170; 0.6666666667] atol = 1e-5

    # ------------------
    # Test model with zero Directions
    # ------------------   
    deaddfmzeros = deaddfm(X, Y, Gx = zeros(size(X)), Gy = zeros(size(Y)), rts = :VRS)
    @test efficiency(deaddfmzeros) == zeros(8)

end
