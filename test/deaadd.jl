# Tests for Radial DEA Models
@testset "AdditiveDEAModel" begin

    ## Test Additive DEA Models with FLS Book data
    # Test against results with R
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    # Additive CRS
    deaaddcrs1 = deaadd(X, Y, :Custom, wX = ones(size(X)), wY = ones(size(Y)), rts = :CRS)

    @test nobs(deaaddcrs1) == 11
    @test ninputs(deaaddcrs1) == 2
    @test noutputs(deaaddcrs1) == 1
    @test efficiency(deaaddcrs1) ≈ [0.0000000000;
                               10.76923077;
                               10.83783784;
                               0.0000000000;
                               22.15384615;
                               17.92307692;
                               0.0000000000;
                               12.07692308;
                               11.61379310;
                               35.000000000;
                               4.000000000]

    @test convert(Matrix, peers(deaaddcrs1)) ≈
    [ 1.0000000000  0  0 0.0000000000  0  0 0.0000000000  0  0   0   0
      0.0000000000  0  0 0.5384615385  0  0 0.0000000000  0  0   0   0
      0.1216216216  0  0 0.9054054054  0  0 0.0000000000  0  0   0   0
      0.0000000000  0  0 1.0000000000  0  0 0.0000000000  0  0   0   0
      0.0000000000  0  0 0.3076923077  0  0 0.0000000000  0  0   0   0
      0.0000000000  0  0 0.3461538462  0  0 0.0000000000  0  0   0   0
      0.0000000000  0  0 0.0000000000  0  0 1.0000000000  0  0   0   0
      0.0000000000  0  0 1.1538461538  0  0 0.0000000000  0  0   0   0
      0.0000000000  0  0 0.4689655172  0  0 0.6965517241  0  0   0   0
      0.0000000000  0  0 1.0000000000  0  0 0.0000000000  0  0   0   0
      1.0000000000  0  0 0.0000000000  0  0 0.0000000000  0  0   0   0]

    # Additive VRS model
    deaaddvrs1 = deaadd(X, Y, :Custom, wX = ones(size(X)), wY = ones(size(Y)), rts = :VRS)

    @test nobs(deaaddvrs1) == 11
    @test ninputs(deaaddvrs1) == 2
    @test noutputs(deaaddvrs1) == 1
    @test efficiency(deaaddvrs1) ≈ [0.0000000000;
                               7.333333333;
                               0.0000000000;
                               0.0000000000;
                               18.000000000;
                               0.0000000000;
                               0.0000000000;
                               0.0000000000;
                               0.0000000000;
                               35.000000000;
                               4.000000000]

    @test convert(Matrix, peers(deaaddvrs1)) ≈
    [1.000000000  0  0  0  0  0 0.00000000000  0  0   0   0;
     0.6666666667 0  0  0  0  0 0.3333333333   0  0   0   0;
     0.000000000  0  1  0  0  0 0.00000000000  0  0   0   0;
     0.000000000  0  0  1  0  0 0.00000000000  0  0   0   0;
     1.000000000  0  0  0  0  0 0.00000000000  0  0   0   0;
     0.000000000  0  0  0  0  1 0.00000000000  0  0   0   0;
     0.000000000  0  0  0  0  0 1.00000000000  0  0   0   0;
     0.000000000  0  0  0  0  0 0.00000000000  1  0   0   0;
     0.000000000  0  0  0  0  0 0.00000000000  0  1   0   0;
     0.000000000  0  0  1  0  0 0.00000000000  0  0   0   0;
     1.000000000  0  0  0  0  0 0.00000000000  0  0   0   0]

    # Check weights equal to one equals default model
    deaadddefaultcrs = deaadd(X, Y, rts = :CRS)
    @test efficiency(deaaddcrs1) ≈ efficiency(deaadddefaultcrs)

    deaadddefaultvrs = deaadd(X, Y)
    @test efficiency(deaaddvrs1) ≈ efficiency(deaadddefaultvrs)

    # Test model :Ones equals model with all weights equal to 1
    deaaddonescrs = deaadd(X, Y, :Ones, rts = :CRS)
    @test deaaddonescrs.weights == :Ones
    @test efficiency(deaaddonescrs) ≈ efficiency(deaadddefaultcrs)

    deaaddonesvrs = deaadd(X, Y, :Ones)
    @test deaaddonesvrs.weights == :Ones
    @test efficiency(deaaddonesvrs) ≈ efficiency(deaadddefaultvrs)

    # MIP CRS
    deaaddmipcrs = deaadd(X, Y, :MIP, rts = :CRS)
    @test deaaddmipcrs.weights == :MIP
    @test efficiency(deaaddmipcrs) ≈ [0.0000000000;
                               0.7577160494;
                               0.4168399168;
                               0.0000000000;
                               2.2219512195;
                               1.1478260870;
                               0.0000000000;
                               0.4867909868;
                               0.4041184041;
                               1.0726153846;
                               0.2352941176]
    @test deaaddmipcrs.slackX ≈ [0 0;
                                 3.037037037 6.814814815;
                                 0 10.837837838;
                                 0 0;
                                 0 0;
                                 8.000000000  0;
                                 0 0;
                                 7.384615385 4.692307692;
                                 8.296296296 2.518518519;
                                 0.000000000  8.200000000;
                                 0 4.000000000]

    @test deaaddmipcrs.slackY ≈ [0;
                                  0;
                                  0;
                                  0;
                                  17.77560976 ;
                                  7.20000000;
                                  0;
                                  0;
                                  0;
                                  19.36000000;
                                  0]

    # MIP VRS
    deaaddmipvrs = deaadd(X, Y, :MIP, rts = :VRS)
    @test deaaddmipvrs.weights == :MIP
    @test efficiency(deaaddmipvrs) ≈ [0.0000000000;
                               0.5075187970;
                               0.0000000000;
                               0.0000000000;
                               2.2039473684;
                               0.0000000000;
                               0.0000000000;
                               0.0000000000;
                               0.0000000000;
                               1.0432234432;
                               0.2352941176]
    @test deaaddmipvrs.slackX ≈ [0 0;
                                 0 0;
                                 0 0;
                                 0 0;
                                 0 0;
                                 0 0;
                                 0 0;
                                 0 0;
                                 0 0;
                                 17 15;
                                 0 4]
    @test deaaddmipvrs.slackY ≈ [0;
                                 7.105263158;
                                 0;
                                 0;
                                 17.631578947;
                                 0;
                                 0;
                                 0;
                                 0;
                                 1.000000000;
                                 0]

    # Normalized CRS
    deaaddnormalizedcrs = deaadd(X, Y, :Normalized, rts = :CRS)
    @test deaaddnormalizedcrs.weights == :Normalized
    @test efficiency(deaaddnormalizedcrs) ≈ [0.0000000000;
                                 1.3569256615;
                                 1.7407259078;
                                 0.0000000000;
                                 2.6877810368;
                                 1.6953153107;
                                 0.0000000000;
                                 1.6540884810;
                                 1.1212042237;
                                 4.0172855145;
                                 0.6424624298]

    @test deaaddnormalizedcrs.slackX ≈ [0.000000000  0.000000000;
                                 3.037037037  6.814814815;
                                 0.000000000 10.837837838;
                                 0.000000000  0.000000000;
                                 10.592592593 11.037037037;
                                 14.666666667  2.666666667;
                                 0.000000000  0.000000000;
                                 0.000000000 10.298429319;
                                 8.296296296  2.518518519;
                                 17.925925926 15.370370370;
                                 0.000000000  4.000000000]

    @test deaaddnormalizedcrs.slackY ≈ [0;
                                  0;
                                  0;
                                  0;
                                  0 ;
                                  0;
                                  0;
                                  0;
                                  0;
                                  0;
                                  0]

    # Normalized VRS
    deaaddnormalizedvrs = deaadd(X, Y, :Normalized, rts = :VRS)
    @test deaaddnormalizedvrs.weights == :Normalized
    @test efficiency(deaaddnormalizedvrs) ≈ [0.0000000000;
                                 0.8049248943;
                                 0.0000000000;
                                 0.0000000000;
                                 2.0149704908;
                                 0.0000000000;
                                 0.0000000000;
                                 0.0000000000;
                                 0.0000000000;
                                 3.9898943952;
                                 0.6424624298]

    @test deaaddnormalizedvrs.slackX ≈ [0  0.00;
                                 0  0.65;
                                 0  0.00;
                                 0  0.00;
                                 0  2.95;
                                 0  0.00;
                                 0  0.00;
                                 0  0.00;
                                 0  0.00;
                                 17 15.00;
                                 0  4.00]

    @test deaaddnormalizedvrs.slackY ≈ [0.00;
                                6.25;
                                0.00;
                                0.00;
                                13.75;
                                0.00;
                                0.00;
                                0.00;
                                0.00;
                                1.00;
                                0.00]

    # RAM CRS
    deaaddramcrs = deaadd(X, Y, :RAM, rts = :CRS)
    @test deaaddramcrs.weights == :RAM
    @test efficiency(deaaddramcrs) ≈ [0.0000000000;
                                0.14094094094;
                                0.18063063063;
                                0.0000000000;
                                0.27937937938;
                                0.17657657658;
                                0.0000000000;
                                0.17164048866;
                                0.11671671672;
                                0.41766766767;
                                0.06666666667]

    @test deaaddramcrs.slackX ≈ [ 0.000000000  0.000000000;
                                3.037037037  6.814814815;
                                0.000000000  10.837837838;
                                0.000000000  0.000000000;
                                10.592592593 11.037037037;
                                14.666666667 2.666666667;
                                0.000000000  0.000000000;
                                0.000000000  10.298429319;
                                8.296296296  2.518518519;
                                17.925925926 15.370370370;
                                0.000000000  4.000000000]

    @test deaaddramcrs.slackY ≈ [0;
                                   0;
                                   0;
                                   0;
                                   0 ;
                                   0;
                                   0;
                                   0;
                                   0;
                                   0;
                                   0]

    # RAM VRS
    deaaddramvrs = deaadd(X, Y, :RAM, rts = :VRS)
    @test deaaddramvrs.weights == :RAM
    @test efficiency(deaaddramvrs) ≈ [0.0000000000;
                              0.10297482838;
                              0.0000000000;
                              0.0000000000;
                              0.25553012967;
                              0.0000000000;
                              0.0000000000;
                              0.0000000000;
                              0.0000000000;
                              0.41764590678 ;
                              0.06666666667]

    @test deaaddramvrs.slackX ≈ [0 0;
                                0 0;
                                0 0;
                                0 0;
                                0 0;
                                0 0;
                                0 0;
                                0 0;
                                0 0;
                                17 15;
                                0 4]
    @test deaaddramvrs.slackY ≈ [0;
                                7.105263158;
                                0;
                                0;
                                17.631578947;
                                0;
                                0;
                                0;
                                0;
                                1.000000000;
                                0]

    # BAM CRS
    deaaddbamcrs = deaadd(X, Y, :BAM, rts = :CRS)
    @test deaaddbamcrs.weights == :BAM
    @test efficiency(deaaddbamcrs) ≈ [0.0000000000;
                                0.4578892372;
                                0.3051750381;
                                0.0000000000;
                                0.6732181854;
                                0.3239568683;
                                0.0000000000;
                                0.5255235602;
                                0.1913580247;
                                0.6902867780;
                                0.1212121212]

    @test deaaddbamcrs.slackX ≈ [0.000000000  0.000000000;
                                4.110344828  6.000000000;
                                0.000000000  0.000000000;
                                0.000000000  0.000000000;
                                13.000000000  8.000000000;
                                17.493670886  0.000000000;
                                0.000000000  0.000000000;
                                0.000000000  9.225130890;
                                8.296296296  2.518518519;
                                13.296296296 13.518518519;
                                0.000000000  4.000000000]
    @test deaaddbamcrs.slackY ≈ [0;
                                0;
                                5.4931506849;
                                0;
                                0.4520547945 ;
                                0;
                                0;
                                1;
                                0;
                                5;
                                0]

    # BAM VRS
    deaaddbamvrs = deaadd(X, Y, :BAM, rts = :VRS)
    @test deaaddbamvrs.weights == :BAM
    @test efficiency(deaaddbamvrs) ≈ [0.0000000000;
                              0.1998936736;
                              0.0000000000;
                              0.0000000000;
                              0.4329710145;
                              0.0000000000;
                              0.0000000000;
                              0.0000000000;
                              0.0000000000;
                              0.5713608345;
                              0.1212121212]

    @test deaaddbamvrs.slackX ≈ [0 0;
                                6.596491228 0;
                                0 0;
                                0 0;
                                13.000000000 1;
                                0 0;
                                0 0;
                                0 0;
                                0 0;
                                5 11;
                                0 4]
    @test deaaddbamvrs.slackY ≈ [0;
                                0;
                                0;
                                0;
                                4;
                                0;
                                0;
                                0;
                                0;
                                5;
                                0]

    # Test model with Custom weights
    deaddcustomcrs = deaadd(X, Y, :Custom, wX = 1 ./ X, wY = 1 ./ Y, rts = :CRS)
    @test deaddcustomcrs.weights == :Custom
    @test efficiency(deaddcustomcrs) ≈ efficiency(deaaddmipcrs)

    deaddcustomvrs = deaadd(X, Y, :Custom, wX = 1 ./ X, wY = 1 ./ Y, rts = :VRS)
    @test deaddcustomvrs.weights == :Custom
    @test efficiency(deaddcustomvrs) ≈ efficiency(deaaddmipvrs)


    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deaadddefaultcrs_ref_eff = zeros(size(X, 1))
    deaadddefaultvrs_ref_eff = zeros(size(X, 1))

    deaaddcustomcrs_ref_eff = zeros(size(X, 1))
    deaaddcustomvrs_ref_eff = zeros(size(X, 1))

    deaaddonescrs_ref_eff = zeros(size(X, 1))
    deaaddonesvrs_ref_eff = zeros(size(X, 1))

    deaaddmipcrs_ref_eff = zeros(size(X, 1))
    deaaddmipvrs_ref_eff = zeros(size(X, 1))

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]

        deaadddefaultcrs_ref_eff[i] = efficiency(deaadd(Xeval, Yeval, rts = :CRS, Xref = Xref, Yref = Yref))[1]
        deaadddefaultvrs_ref_eff[i] = efficiency(deaadd(Xeval, Yeval, rts = :VRS, Xref = Xref, Yref = Yref))[1]

        deaaddcustomcrs_ref_eff[i] = efficiency(deaadd(Xeval, Yeval, :Custom, rts = :CRS, wX = 1 ./ Xeval, wY = 1 ./ Yeval, Xref = Xref, Yref = Yref))[1]
        deaaddcustomvrs_ref_eff[i] = efficiency(deaadd(Xeval, Yeval, :Custom, rts = :VRS, wX = 1 ./ Xeval, wY = 1 ./ Yeval, Xref = Xref, Yref = Yref))[1]

        deaaddonescrs_ref_eff[i] = efficiency(deaadd(Xeval, Yeval, :Ones, rts = :CRS, Xref = Xref, Yref = Yref))[1]
        deaaddonesvrs_ref_eff[i] = efficiency(deaadd(Xeval, Yeval, :Ones, rts = :VRS, Xref = Xref, Yref = Yref))[1]

        deaaddmipcrs_ref_eff[i] = efficiency(deaadd(Xeval, Yeval, :MIP, rts = :CRS, Xref = Xref, Yref = Yref))[1]
        deaaddmipvrs_ref_eff[i] = efficiency(deaadd(Xeval, Yeval, :MIP, rts = :VRS, Xref = Xref, Yref = Yref))[1]

    end

    @test deaadddefaultcrs_ref_eff ≈ efficiency(deaadddefaultcrs)
    @test deaadddefaultvrs_ref_eff ≈ efficiency(deaadddefaultvrs)

    @test deaaddcustomcrs_ref_eff ≈ efficiency(deaaddmipcrs)
    @test deaaddcustomvrs_ref_eff ≈ efficiency(deaaddmipvrs)

    @test deaaddonescrs_ref_eff ≈ efficiency(deaaddcrs1)
    @test deaaddonesvrs_ref_eff ≈ efficiency(deaaddvrs1)

    @test deaaddmipcrs_ref_eff ≈ efficiency(deaaddmipcrs)
    @test deaaddmipvrs_ref_eff ≈ efficiency(deaaddmipvrs)

    # Print
    show(IOBuffer(), deaaddcrs1)

    # Test errors
    @test_throws ErrorException deaadd([1; 2 ; 3], [4 ; 5], :Ones) #  Different number of observations
    @test_throws ErrorException deaadd([1; 2], [4 ; 5], Xref = [1; 2; 3; 4], :Ones) # Different number of observations in reference sets
    @test_throws ErrorException deaadd([1 1; 2 2], [4 4; 5 5], Xref = [1 1 1; 2 2 2], :Ones) # Different number of inputs
    @test_throws ErrorException deaadd([1 1; 2 2], [4 4; 5 5], Yref = [4 4 4; 5 5 5], :Ones) # Different number of inputs
    @test_throws ErrorException deaadd([1; 2; 3], [4; 5; 6], :Ones, rts = :Error) # Invalid returns to scale
    @test_throws ErrorException deaadd([1; 2; 3], [4; 5; 6], :Error) # Invalid model

    @test_throws ErrorException deaadd([1; 2 ; 3], [4 ; 5]) #  Different number of observations
    @test_throws ErrorException deaadd([1; 2], [4 ; 5], Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws ErrorException deaadd([1 1; 2 2], [4 4; 5 5], Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws ErrorException deaadd([1 1; 2 2], [4 4; 5 5], Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ErrorException deaadd([1; 2; 3], [4; 5; 6], rts = :Error) # Invalid returns to scale
    @test_throws ErrorException deaadd([1; 2; 3], [4; 5; 6], :Custom, wX = [1; 2; 3; 4]) # Different size of weights
    @test_throws ErrorException deaadd([1; 2; 3], [4; 5; 6], :Custom, wY = [4; 5; 6; 7]) # Different size of weights
  
end
