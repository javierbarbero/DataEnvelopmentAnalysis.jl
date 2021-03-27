# Tests for Enhanced Russell Graph DEA Models
@testset "EnhancedRussellGraphDEAModel" begin

    # Data for ERG test
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]    

    # Enhanced Russell Graph CRS
    deaergcrs = deaerg(X, Y, rts = :CRS)

    @test typeof(deaergcrs) == EnhancedRussellGraphDEAModel

    @test nobs(deaergcrs) == 8
    @test ninputs(deaergcrs) == 1
    @test noutputs(deaergcrs) == 1
    @test efficiency(deaergcrs) ≈ [0.4; 1.0; 0.8; 0.6; 0.4; 0.4; 0.514285714285650175; 0.2]
    @test efficiency(deaergcrs, :beta) ≈ [0.4; 1.0; 1.0; 0.6; 0.4; 0.4; 0.514285714285650175; 0.2]
    @test convert(Matrix, peers(deaergcrs)) ≈ [
        0 0.5     0    0    0    0    0    0
        0 1.0     0    0    0    0    0    0
        0 1.6     0    0    0    0    0    0
        0 3.0     0    0    0    0    0    0
        0 1.5     0    0    0    0    0    0
        0 3.5     0    0    0    0    0    0
        0 3.5     0    0    0    0    0    0
        0 2.353   0    0    0    0    0    0]

    @test slacks(deaergcrs, :X) ≈ [0; 0; 1.6; 0; 0; 0; 0; 0] 
    @test slacks(deaergcrs, :Y) ≈ [1.5; 0; 0; 6.0; 4.5; 10.5; 8.5; 9.412] 

    @test efficiency(deaerg(targets(deaergcrs, :X), targets(deaergcrs, :Y), rts = :CRS)) ≈ ones(8,1) atol = 1e-6

    # Test default is CRS
    @test efficiency(deaerg(X, Y)) ≈ efficiency(deaergcrs)

    # Enhanced Russell Graph VRS
    deaergvrs = deaerg(X, Y, rts = :VRS)

    @test efficiency(deaergvrs) ≈ [1.0; 1.0; 1.0; 1.0; 0.4; 0.4761904761904763; 0.8571428571428572; 0.2] atol = 1e-10
    @test efficiency(deaergvrs, :beta) ≈ [1.0; 1.0; 1.0; 1.0; 0.6; 1.0; 1.0; 0.4706] 
    @test convert(Matrix, peers(deaergvrs)) ≈ [
        1.0  0     0    0    0    0    0    0
        0    1.0   0    0    0    0    0    0
        0    0     1.0  0    0    0    0    0
        0    0.0   0    1    0    0    0    0
        0    1.0   0    0    0    0    0    0
        0    1/3   2/3  0    0    0    0    0
        0    0     0    1    0    0    0    0
        0    1     0    0    0    0    0    0]

    @test slacks(deaergvrs, :X) ≈ [0; 0; 0; 0; 2.0; 22/3; 2; 5.412] 
    @test slacks(deaergvrs, :Y) ≈ [0; 0; 0; 0; 2.0; 0; 0.0; 2.647] 

    @test efficiency(deaerg(targets(deaergvrs, :X), targets(deaergvrs, :Y), rts = :VRS)) ≈ ones(8,1) atol = 1e-6

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deaerg_crs_ref_eff = zeros(size(X, 1))
    deaerg_vrs_ref_eff = zeros(size(X, 1))

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]

        deaerg_crs_ref_eff[i] = efficiency(deaerg(Xeval, Yeval, rts = :CRS, Xref = Xref, Yref = Yref))[1]
        deaerg_vrs_ref_eff[i] = efficiency(deaerg(Xeval, Yeval, rts = :VRS, Xref = Xref, Yref = Yref))[1]
    end

    @test deaerg_crs_ref_eff ≈ efficiency(deaergcrs)
    @test deaerg_vrs_ref_eff ≈ efficiency(deaergvrs)

    # Print
    show(IOBuffer(), deaergcrs)
    show(IOBuffer(), deaergvrs)

    # ------------------
    # Test errors
    # ------------------
    @test_throws DimensionMismatch deaerg([1; 2 ; 3], [4 ; 5]) #  Different number of observations
    @test_throws DimensionMismatch deaerg([1; 2], [4 ; 5], Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws DimensionMismatch deaerg([1 1; 2 2], [4 4; 5 5], Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws DimensionMismatch deaerg([1 1; 2 2], [4 4; 5 5], Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ArgumentError deaerg([1; 2; 3], [4; 5; 6], rts = :Error) # Invalid returns to scale

    @test_throws ArgumentError efficiency(deaergcrs, :Error) # Invalid efficiency type
   
end
