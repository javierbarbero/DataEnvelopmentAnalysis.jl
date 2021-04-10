# Tests for Reverse DDF DEA Models
@testset "ReverseDDFDEAModel" begin

    # ------------------
    # Input oriented
    # ------------------
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6 8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]

    # Reverse DDF for input :ERG
    russelio = dearussell(X, Y, orient = :Input, rts = :VRS)
    rddfergio = dearddf(X, Y, :ERG, orient = :Input, rts = :VRS)
    
    @test efficiency(rddfergio) ≈ [0.0; 0.0; 0.0; 5/12; 0.6; 1/6; 0.35; 0.4375] atol = 1e-5
    @test targets(russelio, :X) ≈ targets(rddfergio, :X)
    @test targets(russelio, :Y) ≈ targets(rddfergio, :Y)

    # Print
    show(IOBuffer(), rddfergio)

    # ------------------
    # Output oriented
    # ------------------
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]   

    # Reverse DDF for output :ERG
    russeloo = dearussell(X, Y, orient = :Output, rts = :VRS)
    rddfergoo = dearddf(X, Y, :ERG, orient = :Output, rts = :VRS)
    
    @test efficiency(rddfergoo) ≈ [0.0; 0.0; 0.0; 13/28; 4/7; 1/3; 11/35; 37/55] atol = 1e-5
    @test targets(russeloo, :X) ≈ targets(rddfergoo, :X)
    @test targets(russeloo, :Y) ≈ targets(rddfergoo, :Y)

    # Print
    show(IOBuffer(), rddfergoo)

    # ------------------
    # Graph
    # ------------------
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]    

    # Reverse DDF for :ERG :VRS
    erg = deaerg(X, Y, rts = :VRS)
    rddferg = dearddf(X, Y, :ERG, orient = :Graph, rts = :VRS)
    
    @test efficiency(rddferg) ≈ [0.0; 0.0; 0.0; 0.0; 0.6; 11/21; 1/7; 0.8] atol = 1e-5
    @test targets(erg, :X) ≈ targets(rddferg, :X)
    @test targets(erg, :Y) ≈ targets(rddferg, :Y)

    # Reverse DDF for :ERG :VCRS
    ergcrs = deaerg(X, Y, rts = :VRS)
    rddfergcrs = dearddf(X, Y, :ERG, orient = :Graph, rts = :VRS)
    
    @test targets(ergcrs, :X) ≈ targets(rddfergcrs, :X)
    @test targets(ergcrs, :Y) ≈ targets(rddfergcrs, :Y)

    # Default orient is :Graph
    @test efficiency(dearddf(X, Y, :ERG, rts = :VRS)) == efficiency(rddferg)

    # Print
    show(IOBuffer(), rddferg)   

    # Reverse DDF for Modified DDF :VRS
    mddf = deamddf(X, Y, Gx = :Observed, Gy = :Observed, rts = :VRS, slack = false)
    rddfmddf = dearddf(X, Y, :MDDF, Gx = :Observed, Gy = :Observed, rts = :VRS)

    @test efficiency(rddfmddf) ≈ [0.0; 0.0; 0.0; 0.0; 7/6; 4/7; 1/7; 2.549936] atol = 1e-5
    @test targets(rddfmddf, :X) ≈ targets(mddf, :X) atol = 1e-5
    @test targets(rddfmddf, :Y) ≈ targets(mddf, :Y) atol = 1e-5
    
    # Reverse DDF for Modified DDF :CRS
    mddfcrs = deamddf(X, Y, Gx = :Observed, Gy = :Observed, rts = :CRS, slack = false)
    rddfmddfcrs = dearddf(X, Y, :MDDF, Gx = :Observed, Gy = :Observed, rts = :CRS)

    @test targets(rddfmddfcrs, :X) ≈ targets(mddfcrs, :X) atol = 1e-5
    @test targets(rddfmddfcrs, :Y) ≈ targets(mddfcrs, :Y) atol = 1e-5

    # Print
    show(IOBuffer(), rddfmddf)

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    rddf_vrs_ref_eff = zeros(size(X, 1))

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]

        rddf_vrs_ref_eff[i] = efficiency(dearddf(Xeval, Yeval, :ERG; rts = :VRS, Xref = Xref, Yref = Yref))[1]
    end

    @test rddf_vrs_ref_eff ≈ efficiency(rddferg)

    # ------------------
    # Test slacks with different data
    # ------------------   
    X = [1 1; 2 2; 1.5 1]
    Y = [1 1; 2 2; 1 0.5]

    mddfslack = deamddf(X, Y, Gx = :Observed, Gy = :Observed, rts = :VRS)
    rddfslack = dearddf(X, Y, :MDDF, Gx = :Observed, Gy = :Observed, rts = :VRS, slack = true)
    @test slacks(rddfslack, :X) ≈ [0 0; 0 0; 0.5 0] atol = 1e-5
    @test slacks(rddfslack, :Y) ≈ [0 0; 0 0; 0 0.5] atol = 1e-5

    # Print
    show(IOBuffer(), rddfslack)

    # ------------------
    # Test errors
    # ------------------
    @test_throws DimensionMismatch dearddf([1; 2 ; 3], [4 ; 5], :ERG) #  Different number of observations
    @test_throws DimensionMismatch dearddf([1; 2], [4 ; 5], :ERG, Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws DimensionMismatch dearddf([1 1; 2 2], [4 4; 5 5], :ERG, Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws DimensionMismatch dearddf([1 1; 2 2], [4 4; 5 5], :ERG, Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ArgumentError dearddf([1; 2; 3], [4; 5; 6], :Error) # Invalid efficiency measure
    @test_throws ArgumentError dearddf([1; 2; 3], [4; 5; 6], :ERG, rts = :Error) # Invalid returns to scale
    @test_throws ArgumentError dearddf([1; 2; 3], [4; 5; 6], :ERG, orient = :Error) # Invalid orientation
    @test_throws ArgumentError dearddf([1; 2; 3], [4; 5; 6], :ERG, Gx = :Ones) # No direction in :ERG
    @test_throws ArgumentError dearddf([1; 2; 3], [4; 5; 6], :MDDF, orient = :Input, Gx = :Ones, Gy = :Ones) # Invalid direction in :MDDF
    @test_throws DimensionMismatch dearddf([1 1; 2 2; 3 3], [4; 5; 6], :MDDF, Gx = [1 1 1; 2 2 2; 3 3 3], Gy = [4; 5; 6]) # Different size of inputs direction
    @test_throws DimensionMismatch dearddf([1; 2; 3], [4 4; 5 5; 6 6], :MDDF, Gx = [1; 2; 3], Gy = [4 4 4; 5 5 5; 6 6 6]) # Different size of inputs direction
    @test_throws ArgumentError dearddf([1; 2; 3], [1; 2; 3], :MDDF, Gx = :Error, Gy = :Ones) # Invalid inuts direction
    @test_throws ArgumentError dearddf([1; 2; 3], [1; 2; 3], :MDDF, Gx = :Ones, Gy = :Error) # Invalid outputs direction

end
