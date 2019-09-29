# Tests for Generalized DF DEA Models
@testset "GeneralizedDFDEAModel" begin

    ## Test Generalized DF DEA Model with Zofío and Prieto (2006) data
    X = [5 3; 2 4; 4 2; 4 8; 7 9]
    Y = [7 4; 10 8; 8 10; 5 4; 3 6]

    # alpha = 0.5 CRS equals Input Oriented CRS
    deaio = dea(X, Y, orient = :Input, rts = :CRS)
    deagdf05crs = deagdf(X, Y, 0.5, rts = :CRS)
    @test efficiency(deaio) ≈ efficiency(deagdf05crs) atol = 1e-7

    @test nobs(deagdf05crs) == 5
    @test ninputs(deagdf05crs) == 2
    @test noutputs(deagdf05crs) == 2

    # alphpa = 0.5 VRS
    deagdf05vrs = deagdf(X, Y, 0.5, rts = :VRS)

    @test nobs(deagdf05vrs) == 5
    @test ninputs(deagdf05vrs) == 2
    @test noutputs(deagdf05vrs) == 2
    @test efficiency(deagdf05vrs) ≈ [0.682;
                               1;
                               1;
                               0.250;
                               0.360] atol = 1e-3

    # Test no slacks
    deagdfnoslack = deagdf(X, Y, 0.5, rts = :VRS, slack = false)
    @test efficiency(deagdfnoslack) == efficiency(deagdf05vrs)
    @test isempty(slacks(deagdfnoslack, :X)) == 1
    @test isempty(slacks(deagdfnoslack, :Y)) == 1

     ## Test if one-by-one DEA using evaluation and reference sets match initial results
     deagdf05crs_ref_eff = zeros(size(X, 1))

     deagdf05vrs_ref_eff = zeros(size(X, 1))

     Xref = X[:,:]
     Yref = Y[:,:]

     for i = 1:size(X, 1)
         Xeval = X[i:i,:]
         Xeval = Xeval[:,:]
         Yeval = Y[i:i,:]
         Yeval = Yeval[:,:]

         deagdf05crs_ref_eff[i] = efficiency(deagdf(Xeval, Yeval, 0.5, rts = :CRS, Xref = Xref, Yref = Yref))[1]

         deagdf05vrs_ref_eff[i] = efficiency(deagdf(Xeval, Yeval, 0.5, rts = :VRS, Xref = Xref, Yref = Yref))[1]
     end

     @test deagdf05crs_ref_eff ≈ efficiency(deagdf05crs)
     @test deagdf05vrs_ref_eff ≈ efficiency(deagdf05vrs)

    # Print
    show(IOBuffer(), deagdf05crs)
    show(IOBuffer(), deagdfnoslack)

    # Test errors
    @test_throws ErrorException deagdf([1; 2 ; 3], [4 ; 5], 0.5) #  Different number of observations
    @test_throws ErrorException deagdf([1; 2], [4 ; 5], 0.5, Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws ErrorException deagdf([1 1; 2 2], [4 4; 5 5], 0.5, Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws ErrorException deagdf([1 1; 2 2], [4 4; 5 5], 0.5, Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ErrorException deagdf([1; 2; 3], [4; 5; 6], 0.5, rts = :Error) # Invalid returns to scale


end
