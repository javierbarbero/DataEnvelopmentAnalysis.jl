# Tests for Radial Multiplier DEA Models
@testset "RadialMultiplierDEAModel" begin

    ## Test Radial Multiplier DEA Models with FLS Book data
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    # Input oriented CRS
    deamio = deam(X, Y, orient = :Input, rts = :CRS)

    @test typeof(deamio) == RadialMultiplierDEAModel

    @test nobs(deamio) == 11
    @test ninputs(deamio) == 2
    @test noutputs(deamio) == 1
    @test efficiency(deamio) ≈ [
        1.0000000000;
        0.6222896791;
        0.8198562444;
        1.0000000000;
        0.3103709311;
        0.5555555556;
        1.0000000000;
        0.7576690896;
        0.8201058201;
        0.4905660377;
        1.0000000000] atol = 1e-5

    @test multipliers(deamio, :X) ≈ [
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
    @test multipliers(deamio, :Y) ≈ [
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

    @test efficiency(deam(targets(deamio, :X), targets(deamio, :Y), orient = :Input, rts = :CRS)) ≈ ones(11) atol = 1e-5
    @test rts(deamio) ≈ zeros(11)

    # Otuput oriented CRS
    deamoo = deam(X, Y, orient = :Output, rts = :CRS)

    @test nobs(deamoo) == 11
    @test ninputs(deamoo) == 2
    @test noutputs(deamoo) == 1
    @test efficiency(deamoo) ≈ [
        1.0000000000;
        1.606968641;
        1.219726027;
        1.0000000000;
        3.221951220;
        1.800000000;
        1.0000000000;
        1.319837398;
        1.219354839;
        2.038461538;
        1.0000000000] atol = 1e-5

    @test multipliers(deamoo, :X) ≈ [
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
    @test multipliers(deamoo, :Y) ≈ [
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

    @test efficiency(deam(targets(deamoo, :X), targets(deamoo, :Y), orient = :Output, rts = :CRS)) ≈ ones(11) atol = 1e-5
    @test rts(deamoo) ≈ zeros(11)

    # Input oriented VRS
    deamiovrs = deam(X, Y, orient = :Input, rts = :VRS)

    @test nobs(deamiovrs) == 11
    @test ninputs(deamiovrs) == 2
    @test noutputs(deamiovrs) == 1
    @test efficiency(deamiovrs) ≈ [
        1.0000000000;
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

    @test multipliers(deamiovrs, :X) ≈ [
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
    @test multipliers(deamiovrs, :Y) ≈ [
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

    @test efficiency(deam(targets(deamiovrs, :X), targets(deamiovrs, :Y), orient = :Input, rts = :VRS)) ≈ ones(11) atol = 1e-5
    @test rts(deamiovrs) ≈ [0.0; -0.661826; 0.314757; 0.0; -0.711640; -0.666667; 0.0; 1.626970; 1.214286; -0.135409; 0.0] atol = 1e-5

    # Output oriented VRS
    deamoovrs = deam(X, Y, orient = :Output, rts = :VRS)

    @test nobs(deamoovrs) == 11
    @test ninputs(deamoovrs) == 2
    @test noutputs(deamoovrs) == 1
    @test efficiency(deamoovrs) ≈ [
        1.0000000000;
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

    @test multipliers(deamoovrs, :X) ≈ [
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
    @test multipliers(deamoovrs, :Y) ≈ [
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

    @test efficiency(deam(targets(deamoovrs, :X), targets(deamoovrs, :Y), orient = :Output, rts = :VRS)) ≈ ones(11) atol = 1e-5
    @test rts(deamoovrs) ≈ [0.0; -0.703008; 0.239403; 0.0; -1.230263; -2.0; 0.0; 0.619333; 0.548387; 1.192308; 0.0] atol = 1e-5

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deamio_ref_eff = zeros(size(X, 1))
    deamoo_ref_eff = zeros(size(X, 1))

    deamiovrs_ref_eff = zeros(size(X, 1))
    deamoovrs_ref_eff = zeros(size(X, 1))

    deamiovrs_ref_multX = zeros(size(X))
    deamiovrs_ref_multY = zeros(size(Y))

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]

        deamio_ref_eff[i] = efficiency(deam(Xeval, Yeval, orient = :Input, rts = :CRS, Xref = Xref, Yref = Yref))[1]
        deamoo_ref_eff[i] = efficiency(deam(Xeval, Yeval, orient = :Output, rts = :CRS, Xref = Xref, Yref = Yref))[1]

        deamiovrs_ref_eff[i] = efficiency(deam(Xeval, Yeval, orient = :Input, rts = :VRS, Xref = Xref, Yref = Yref))[1]
        deamoovrs_ref_eff[i] = efficiency(deam(Xeval, Yeval, orient = :Output, rts = :VRS, Xref = Xref, Yref = Yref))[1]

        deamiovrs_ref_multX[i,:] = multipliers(deam(Xeval, Yeval, orient = :Input, rts = :VRS, Xref = Xref, Yref = Yref), :X)
        deamiovrs_ref_multY[i,:] = multipliers(deam(Xeval, Yeval, orient = :Input, rts = :VRS, Xref = Xref, Yref = Yref), :Y)
    end

    @test deamio_ref_eff ≈ efficiency(deamio)
    @test deamoo_ref_eff ≈ efficiency(deamoo)

    @test deamiovrs_ref_eff ≈ efficiency(deamiovrs)
    @test deamoovrs_ref_eff ≈ efficiency(deamoovrs)

    @test deamiovrs_ref_multX ≈ multipliers(deamiovrs, :X) 
    @test deamiovrs_ref_multY ≈ multipliers(deamiovrs, :Y) 

    # Print
    show(IOBuffer(), deamio)

    # Test errors
    @test_throws DimensionMismatch deam([1; 2 ; 3], [4 ; 5]) #  Different number of observations
    @test_throws DimensionMismatch deam([1; 2], [4 ; 5], Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws DimensionMismatch deam([1 1; 2 2], [4 4; 5 5], Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws DimensionMismatch deam([1 1; 2 2], [4 4; 5 5], Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ArgumentError deam([1; 2; 3], [4; 5; 6], orient = :Error) # Invalid orientation
    @test_throws ArgumentError deam([1; 2; 3], [4; 5; 6], rts = :Error) # Invalid returns to scale

    @test_throws ArgumentError targets(deamio, :Error)    # Invalid target
    @test_throws ArgumentError multipliers(deamio, :Error) # Invalid slacks

    # ------------------
    # Test Vector and Matrix inputs and outputs
    # ------------------
    # Tests against results in R

    # Inputs is Matrix, Outputs is Vector
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deam(X, Y, orient = :Input)) ≈ [1; 1; 1; 0.6; 0.4; 1; 0.6666666667; 0.625]

    # Inputs is Vector, Output is Matrix
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]

    @test efficiency(deam(X, Y, orient = :Output)) ≈ [1; 1; 1; 1.555555556; 2.333333333; 1; 1.272727273; 1.6]

    # Inputs is Vector, Output is Vector
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]

    @test efficiency(deam(X, Y, orient = :Input)) ≈ [0.4; 1; 0.8; 0.6; 0.4; 0.4; 0.5142857143; 0.2]

end
