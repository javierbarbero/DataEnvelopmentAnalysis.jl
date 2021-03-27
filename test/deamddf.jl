# Tests for Modified DDF DEA Models
@testset "ModifiedDDFDEAModel" begin

    # Data for Modified DDF test
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]    

    # Modified DDF CRS Ones
    deamddfcrs = deamddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :CRS)

    @test efficiency(deamddfcrs) ≈ [1.5; 0; 2; 6; 4.5; 10.5; 8.5; 9.412] atol = 1e-5
    @test efficiency(deamddfcrs, :X) ≈ zeros(8, 1) atol = 1e-5
    @test efficiency(deamddfcrs, :Y) ≈ [1.5; 0; 2; 6; 4.5; 10.5; 8.5; 9.412] atol = 1e-5
    @test convert(Matrix, peers(deamddfcrs)) ≈ [
        0 0.5     0    0    0    0    0    0
        0 1.0     0    0    0    0    0    0
        0 2.0     0    0    0    0    0    0
        0 3.0     0    0    0    0    0    0
        0 1.5     0    0    0    0    0    0
        0 3.5     0    0    0    0    0    0
        0 3.5     0    0    0    0    0    0
        0 2.353   0    0    0    0    0    0]  atol = 1e-5

    @test slacks(deamddfcrs, :X) ≈ zeros(8, 1)
    @test slacks(deamddfcrs, :Y) ≈ zeros(8, 1)

    @test efficiency(deamddf(targets(deamddfcrs, :X), targets(deamddfcrs, :Y), Gx = :Ones, Gy = :Ones, rts = :CRS)) ≈ zeros(8,1) atol = 1e-5

    # Test no slacks
    deamddfnoslack = deamddf(X, Y, Gx = :Ones, Gy = :Ones, slack = false)
    @test efficiency(deamddfnoslack) == efficiency(deamddfcrs)
    @test isempty(slacks(deamddfnoslack, :X)) == 1
    @test isempty(slacks(deamddfnoslack, :Y)) == 1

    # Test default is CRS and :Ones
    @test efficiency(deamddf(X, Y, Gx = :Ones, Gy = :Ones)) ≈ efficiency(deamddfcrs)

    # Modified DDF VRS Ones
    deamddfvrs = deamddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :VRS)

    @test efficiency(deamddfvrs) ≈ [0; 0; 0; 0; 4; 22/3; 2; 8.059] atol = 1e-5
    @test efficiency(deamddfvrs, :X) ≈ [0; 0; 0; 0; 2.0; 22/3; 2; 5.412] atol = 1e-5 
    @test efficiency(deamddfvrs, :Y) ≈ [0; 0; 0; 0; 2; 0; 0; 2.647] atol = 1e-5 
    @test convert(Matrix, peers(deamddfvrs)) ≈ [
        1.0  0     0    0    0    0    0    0
        0    1.0   0    0    0    0    0    0
        0    0     1.0  0    0    0    0    0
        0    0.0   0    1    0    0    0    0
        0    1.0   0    0    0    0    0    0
        0    1/3   2/3  0    0    0    0    0
        0    0     0    1    0    0    0    0
        0    1     0    0    0    0    0    0] atol = 1e-5 

    @test slacks(deamddfvrs, :X) ≈ zeros(8, 1)
    @test slacks(deamddfvrs, :Y) ≈  zeros(8, 1)

    @test efficiency(deamddf(targets(deamddfvrs, :X), targets(deamddfvrs, :Y), Gx = :Ones, Gy = :Ones, rts = :VRS)) ≈ zeros(8,1) atol = 1e-5

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deamddf_crs_ref_eff = zeros(size(X, 1))
    deamddf_vrs_ref_eff = zeros(size(X, 1))

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]

        deamddf_crs_ref_eff[i] = efficiency(deamddf(Xeval, Yeval, Gx = :Ones, Gy = :Ones,  rts = :CRS, Xref = Xref, Yref = Yref))[1]
        deamddf_vrs_ref_eff[i] = efficiency(deamddf(Xeval, Yeval, Gx = :Ones, Gy = :Ones, rts = :VRS, Xref = Xref, Yref = Yref))[1]
    end

    @test deamddf_crs_ref_eff ≈ efficiency(deamddfcrs)
    @test deamddf_vrs_ref_eff ≈ efficiency(deamddfvrs)

    # Modified DDF VRS Observed
    deamddfvrsobs = deamddf(X, Y, Gx = :Observed, Gy = :Observed, rts = :VRS)

    @test efficiency(deamddfvrsobs) ≈ [0; 0; 0; 0; 35/30; 0.571429; 0.142857; 2.549936] atol = 1e-5
    @test efficiency(deamddfvrsobs, :X) ≈ [0; 0; 0; 0; 0; 0.428571; 0.142857; 0.108763] atol = 1e-5
    @test efficiency(deamddfvrsobs, :Y) ≈ [0; 0; 0; 0; 35/30; 0.142857; 0.0; 2.441174] atol = 1e-5

    # Modified DDF VRS Mean
    deamddfvrsmean = deamddf(X, Y, Gx = :Mean, Gy = :Mean, rts = :VRS)

    @test efficiency(deamddfvrsmean) ≈ [0; 0; 0; 0; 0.631299; 0.871894; 0.230508; 1.181294] atol = 1e-5
    @test efficiency(deamddfvrsmean, :X) ≈ [0; 0; 0; 0; 0; 0.691523; 0.230508; 0.162738] atol = 1e-5
    @test efficiency(deamddfvrsmean, :Y) ≈ [0; 0; 0; 0; 0.631299; 0.180371; 0.0; 1.018556] atol = 1e-5

    # Test Custom direction
    @test efficiency(deamddf(X, Y, Gx = ones(8,1), Gy = ones(8,1), rts = :VRS)) ≈ efficiency(deamddfvrs)

    # Print
    show(IOBuffer(), deamddfcrs)
    show(IOBuffer(), deamddfvrs)

    # ------------------
    # Test slacks with different data
    # ------------------   
    X = [1 1; 2 2; 1.5 1]
    Y = [1 1; 2 2; 1 0.5]

    deamddfslack = deamddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :VRS)
    @test slacks(deamddfslack, :X) ≈ [0 0; 0 0; 0.5 0] atol = 1e-5
    @test slacks(deamddfslack, :Y) ≈ [0 0; 0 0; 0 0.5] atol = 1e-5

    # ------------------
    # Test errors
    # ------------------
    @test_throws DimensionMismatch deamddf([1; 2 ; 3], [4 ; 5], Gx = :Ones, Gy = :Ones) #  Different number of observations
    @test_throws DimensionMismatch deamddf([1; 2], [4 ; 5], Gx = :Ones, Gy = :Ones, Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws DimensionMismatch deamddf([1 1; 2 2], [4 4; 5 5], Gx = :Ones, Gy = :Ones, Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws DimensionMismatch deamddf([1 1; 2 2], [4 4; 5 5], Gx = :Ones, Gy = :Ones, Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ArgumentError deamddf([1; 2; 3], [4; 5; 6], Gx = :Ones, Gy = :Ones, rts = :Error) # Invalid returns to scale
    @test_throws DimensionMismatch deamddf([1 1; 2 2; 3 3], [4; 5; 6], Gx = [1 1 1; 2 2 2; 3 3 3], Gy = [4; 5; 6]) # Different size of inputs direction
    @test_throws DimensionMismatch deamddf([1; 2; 3], [4 4; 5 5; 6 6], Gx = [1; 2; 3], Gy = [4 4 4; 5 5 5; 6 6 6]) # Different size of inputs direction
    @test_throws ArgumentError deamddf([1; 2; 3], [1; 2; 3], Gx = :Error, Gy = :Ones) # Invalid inuts direction
    @test_throws ArgumentError deamddf([1; 2; 3], [1; 2; 3], Gx = :Ones, Gy = :Error) # Invalid outputs direction

    @test_throws ArgumentError efficiency(deamddfvrs, :Error) # Invalid efficiency type
   
end
