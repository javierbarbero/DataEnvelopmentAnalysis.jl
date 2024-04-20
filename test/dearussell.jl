# Tests for Russell DEA Models
@testset "RussellDEAModel" begin

    # ------------------
    # Input oriented
    # ------------------
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6 8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]

    # Input oriented CRS
    dearussellio = dearussell(X, Y, orient = :Input, rts = :CRS)

    @test typeof(dearussellio) == RussellDEAModel

    @test nobs(dearussellio) == 8
    @test ninputs(dearussellio) == 2
    @test noutputs(dearussellio) == 1
    @test efficiency(dearussellio, :X) ≈ [1 1; 1 1; 1 1; 0.5 2/3; 0.4 0.4; 2/3 1; 0.5 0.8; 0.625 0.5]
    @test_throws ArgumentError efficiency(dearussellio, :Y) # No output efficiency in input oriented model 
    @test efficiency(dearussellio) ≈ [1; 1; 1; 0.5833333333333334; 0.4; 0.8333333333333334; 0.65; 0.5625]
    @test convert(Matrix, peers(dearussellio)) == 
            [ 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0]
    @test slacks(dearussellio, :X) ≈ zeros(8,2) atol = 1e-10
    @test slacks(dearussellio, :Y) ≈ zeros(8,1) atol = 1e-10

    @test efficiency(dearussell(targets(dearussellio, :X), targets(dearussellio, :Y), orient = :Input, rts = :CRS)) ≈ ones(8,1)

    # Test default is :Input with :CRS
    efficiency(dearussell(X, Y)) == efficiency(dearussellio)

    # Test no slacks
    dearussellionoslack = dearussell(X, Y, slack = false)
    @test efficiency(dearussellionoslack) == efficiency(dearussellio)
    @test isempty(slacks(dearussellionoslack, :X)) == 1
    @test isempty(slacks(dearussellionoslack, :Y)) == 1

    # Input oriented VRS
    dearusselliovrs  = dearussell(X, Y, orient = :Input, rts = :VRS)
    @test efficiency(dearusselliovrs) ≈ [1; 1; 1; 0.5833333333333334; 0.4; 0.8333333333333334; 0.65; 0.5625]

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    dearussell_io_ref_eff = zeros(size(X, 1))
    dearussell_iovrs_ref_eff = zeros(size(X, 1))

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]

        dearussell_io_ref_eff[i] = efficiency(dearussell(Xeval, Yeval, orient = :Input, rts = :CRS, Xref = Xref, Yref = Yref))[1]
        dearussell_iovrs_ref_eff[i] = efficiency(dearussell(Xeval, Yeval, orient = :Input, rts = :VRS, Xref = Xref, Yref = Yref))[1]
    end

    @test dearussell_io_ref_eff == efficiency(dearussellio)
    @test dearussell_iovrs_ref_eff == efficiency(dearusselliovrs)

    # Add new DMU with slack 
    dearussellioslack = dearussell([X; 2 2], [Y; 2], orient = :Input, rts = :VRS)
    @test slacks(dearussellioslack, :Y) ≈ [1; 0; 0; 1; 1; 0; 0; 0; 0]

    # Print
    show(IOBuffer(), dearussellio)
    show(IOBuffer(), dearussellionoslack)
    
    # ------------------
    # Output oriented
    # ------------------
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]    

    # Output oriented CRS
    dearusselloo = dearussell(X, Y, orient = :Output, rts = :CRS)

    @test nobs(dearusselloo) == 8
    @test ninputs(dearusselloo) == 1
    @test noutputs(dearusselloo) == 2
    @test efficiency(dearusselloo, :Y) ≈ [1 1; 1 1; 1 1; 7/3 1.4; 7/3 7/3; 1 2; 7/6 7/4; 46/9 1]
    @test_throws ArgumentError efficiency(dearusselloo, :X) # No input efficiency in input oriented model 
    @test efficiency(dearusselloo) ≈ [1.0; 1.0; 1.0; 1.8666666666666665; 2.333333333333333; 1.5; 1.4583333333333335; 3.0555555555555554]
    @test convert(Matrix, peers(dearusselloo)) ≈ 
            [ 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            1/3  0.0  2/3  0.0  0.0  0.0  0.0  0.0] atol = 1e-10
    @test slacks(dearusselloo, :X) ≈ zeros(8,1)
    @test slacks(dearusselloo, :Y) ≈ zeros(8,2) atol = 1e-14

    @test efficiency(dearussell(targets(dearusselloo, :X), targets(dearusselloo, :Y), orient = :Output, rts = :CRS)) ≈ ones(8,1)

    # Test no slacks
    dearusselloonoslack = dearussell(X, Y, orient = :Output, slack = false)
    @test efficiency(dearusselloonoslack) == efficiency(dearusselloo)
    @test isempty(slacks(dearusselloonoslack, :X)) == 1
    @test isempty(slacks(dearusselloonoslack, :Y)) == 1

    # Output oriented VRS
    dearusselloovrs  = dearussell(X, Y, orient = :Output, rts = :VRS)
    @test efficiency(dearusselloovrs) ≈ [1.0; 1.0; 1.0; 1.8666666666666665; 2.333333333333333; 1.5; 1.4583333333333335; 3.0555555555555554]

    # Output oriented FDH
    dearussellfdh = dearussell(X, Y, orient = :Output, rts = :FDH)

    @test efficiency(dearussellfdh, :Y) ≈ [1 1; 1 1; 1 1; 7/3 1.4; 7/3 7/3; 1 2; 7/6 7/4; 14/3 1.4]
    @test efficiency(dearussellfdh) ≈ [1.0; 1.0; 1.0; 1.8666666666666665; 2.333333333333333; 1.5; 1.4583333333333335; 3.033333333333333]
    @test convert(Matrix, peers(dearussellfdh)) ≈ 
            [ 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
            1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0] atol = 1e-10
    @test slacks(dearussellfdh, :X) ≈ zeros(8,1)
    @test slacks(dearussellfdh, :Y) ≈ zeros(8,2) atol = 1e-14

    @test efficiency(dearussell(targets(dearussellfdh, :X), targets(dearussellfdh, :Y), orient = :Output, rts = :FDH)) ≈ ones(8,1)

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    dearussell_oo_ref_eff = zeros(size(X, 1))
    dearussell_oovrs_ref_eff = zeros(size(X, 1))

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]

        dearussell_oo_ref_eff[i] = efficiency(dearussell(Xeval, Yeval, orient = :Output, rts = :CRS, Xref = Xref, Yref = Yref))[1]
        dearussell_oovrs_ref_eff[i] = efficiency(dearussell(Xeval, Yeval, orient = :Output, rts = :VRS, Xref = Xref, Yref = Yref))[1]
    end

    @test dearussell_oo_ref_eff == efficiency(dearusselloo)
    @test dearussell_oovrs_ref_eff == efficiency(dearusselloovrs)

    # Add new DMU with slack 
    dearussellooslack = dearussell([X; 0.5], [Y; 7 7], orient = :Output, rts = :VRS)
    @test slacks(dearussellooslack, :X) ≈ [0.5; 0; 0; 0.5; 0.5; 0; 0.5; 1/6; 0]

    # Print
    show(IOBuffer(), dearusselloo)
    show(IOBuffer(), dearusselloonoslack)

    # ------------------
    # Graph
    # ------------------
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]    

    # Graph CRS
    dearussellgr = dearussell(X, Y, orient = :Graph, rts = :CRS)

    @test nobs(dearussellgr) == 8
    @test ninputs(dearussellgr) == 1
    @test noutputs(dearussellgr) == 1
    @test efficiency(dearussellgr, :X) ≈ [ 0.632455;
        1.0;
        0.894427;
        0.774597;
        0.632456;
        0.632456;
        0.717137;
        0.447214] atol = 1e-5
    @test efficiency(dearussellgr, :Y) ≈ [ 1.581139;
        1.0;
        1.118034;
        1.290994;
        1.581139;
        1.581139;
        1.394433;
        2.236068] atol = 1e-5
    @test efficiency(dearussellgr) ≈ [ 0.632455;
        1.0;
        0.894427;
        0.774597;
        0.632455;
        0.632455;
        0.717137;
        0.447214] atol = 1e-5

    @test isempty(slacks(dearussellgr, :X)) == 1
    @test isempty(slacks(dearussellgr, :Y)) == 1

    @test efficiency(dearussell(targets(dearussellgr, :X), targets(dearussellgr, :Y), orient = :Graph, rts = :CRS)) ≈ ones(8,1) atol = 1e-6

    # Graph VRS
    dearussellgrvrs  = dearussell(X, Y, orient = :Graph, rts = :VRS)

    @test efficiency(dearussellgrvrs, :X) ≈ [  1.0;
        1.0;
        1.0;
        1.0;
        2/3;
        0.571429;
        0.857143;
        0.424989] atol = 1e-5
    @test efficiency(dearussellgrvrs, :Y) ≈ [ 1.0;
        1.0;
        1.0;
        1.0;
        5/3;
        1.142857;
        1.0;
        2.124947] atol = 1e-5
    @test efficiency(dearussellgrvrs) ≈ [ 1.0;
        1.0;
        1.0;
        1.0;
        0.633333;
        0.723214;
        0.928571;
        0.447795] atol = 1e-5

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    dearussell_gr_ref_eff = zeros(size(X, 1))
    dearussell_grvrs_ref_eff = zeros(size(X, 1))

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]

        dearussell_gr_ref_eff[i] = efficiency(dearussell(Xeval, Yeval, orient = :Graph, rts = :CRS, Xref = Xref, Yref = Yref))[1]
        dearussell_grvrs_ref_eff[i] = efficiency(dearussell(Xeval, Yeval, orient = :Graph, rts = :VRS, Xref = Xref, Yref = Yref))[1]
    end

    @test dearussell_gr_ref_eff == efficiency(dearussellgr)
    @test dearussell_grvrs_ref_eff == efficiency(dearussellgrvrs)

    # Print
    show(IOBuffer(), dearussellgr)

    # ------------------
    # Test errors
    # ------------------
    @test_throws DimensionMismatch dearussell([1; 2 ; 3], [4 ; 5]) #  Different number of observations
    @test_throws DimensionMismatch dearussell([1; 2], [4 ; 5], Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws DimensionMismatch dearussell([1 1; 2 2], [4 4; 5 5], Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws DimensionMismatch dearussell([1 1; 2 2], [4 4; 5 5], Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ArgumentError dearussell([1; 2; 3], [4; 5; 6], orient = :Error) # Invalid orientation
    @test_throws ArgumentError dearussell([1; 2; 3], [4; 5; 6], rts = :Error) # Invalid returns to scale
    @test_throws ArgumentError dearussell([1; 2; 3], [4; 5; 6], orient = :Graph, rts = :FDH) # FDH not implemented for orientation :Graph

    @test_throws ArgumentError efficiency(dearussellio, :Error) # Invalid efficiency type
    

end
