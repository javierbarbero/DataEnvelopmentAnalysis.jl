# Tests for Profitability DEA Models
@testset "ProfitabilityDEAModel" begin

    ## Test Profitability DEA Model with Zofío and Prieto (2006) data
    X = [5 3; 2 4; 4 2; 4 8; 7 9]
    Y = [7 4; 10 8; 8 10; 5 4; 3 6]
    W = [2 1; 2 1; 2 1; 2 1; 2 1.0]
    P = [3 2; 3 2; 3 2; 3 2; 3 2.0]

    # alpha = 0.5 CRS equals Input Oriented CRS
    deaprofbl = deaprofitability(X, Y, W, P)
    @test efficiency(deaprofbl) ≈ [0.388;
                                   1.000;
                                   0.765;
                                   0.250;
                                   0.159] atol = 1e-3
    @test efficiency(deaprofbl, :CRS) ≈ [0.636;
                                1.000;
                                1.000;
                                0.250;
                                0.261] atol = 1e-3
    @test efficiency(deaprofbl, :VRS) ≈ [0.682;
                                1.000;
                                1.000;
                                0.250;
                                0.360] atol = 1e-3
    @test efficiency(deaprofbl, :Scale) ≈ [0.933;
                                1.000;
                                1.000;
                                1.000;
                                0.725] atol = 1e-3
    @test efficiency(deaprofbl, :Allocative) ≈ [0.610;
                                1.000;
                                0.765;
                                1.000;
                                0.609] atol = 1e-3

    # Check defaults
    @test efficiency(deaprofitability(X, Y, W, P, alpha = 0.5)) == efficiency(deaprofbl)
    @test efficiency(deaprofbl, :Economic) == efficiency(deaprofbl)

     ## Test if one-by-one DEA using evaluation and reference sets match initial results
     deaprofbl_ref_eff = zeros(size(X, 1))
     deaprofbl_ref_crs = zeros(size(X, 1))
     deaprofbl_ref_vrs = zeros(size(X, 1))
     deaprofbl_ref_scale = zeros(size(X, 1))
     deaprofbl_ref_alloc = zeros(size(X, 1))
     Xref = X[:,:]
     Yref = Y[:,:]
     Wref = W[:,:]
     Pref = P[:,:]

     for i = 1:size(X, 1)
         Xeval = X[i:i,:]
         Xeval = Xeval[:,:]
         Yeval = Y[i:i,:]
         Yeval = Yeval[:,:]
         Weval = W[i:i,:]
         Weval = Weval[:,:]
         Peval = P[i:i,:]
         Peval = Peval[:,:]

         deaprofbl_ref_eff[i] = efficiency(deaprofitability(Xeval, Yeval, Weval, Peval,  Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref))[1]
         deaprofbl_ref_crs[i] = efficiency(deaprofitability(Xeval, Yeval, Weval, Peval,  Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref), :CRS)[1]
         deaprofbl_ref_vrs[i] = efficiency(deaprofitability(Xeval, Yeval, Weval, Peval,  Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref), :VRS)[1]
         deaprofbl_ref_scale[i] = efficiency(deaprofitability(Xeval, Yeval, Weval, Peval,  Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref), :Scale)[1]
         deaprofbl_ref_alloc[i] = efficiency(deaprofitability(Xeval, Yeval, Weval, Peval,  Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref), :Allocative)[1]

     end

     @test deaprofbl_ref_eff ≈ efficiency(deaprofbl)
     @test deaprofbl_ref_crs ≈ efficiency(deaprofbl, :CRS)
     @test deaprofbl_ref_vrs ≈ efficiency(deaprofbl, :VRS)
     @test deaprofbl_ref_scale ≈ efficiency(deaprofbl, :Scale)
     @test deaprofbl_ref_alloc ≈ efficiency(deaprofbl, :Allocative)

    # Print
    show(IOBuffer(), deaprofbl)

    # Test errors
    @test_throws ErrorException deaprofitability([1; 2 ; 3], [4 ; 5], [1; 1; 1], [4; 5]) #  Different number of observations
    @test_throws ErrorException deaprofitability([1; 2], [4 ; 5], [1; 1], [4 ; 5], Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws ErrorException deaprofitability([1 1; 2 2], [4 4; 5 5], [1 1; 2 2], [4 4; 5 5], Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws ErrorException deaprofitability([1 1; 2 2], [4 4; 5 5], [1 1; 2 2], [4 4; 5 5], Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ErrorException deaprofitability([1; 2; 3], [4; 5; 6], [1; 2; 3; 4], [4; 5; 6]) # Different number of observation in input prices
    @test_throws ErrorException deaprofitability([1; 2; 3], [4; 5; 6], [1; 2; 3], [4; 5; 6; 7]) # Different number of observation in output prices
    @test_throws ErrorException deaprofitability([1; 2; 3], [4 4; 5 5; 6 6], [1; 2; 3], [4 4 4; 5 5 5; 6 6 6]) # Different number of input prices and inputs
    @test_throws ErrorException deaprofitability([1; 2; 3], [4 4; 5 5; 6 6], [1; 2; 3], [4 4 4; 5 5 5; 6 6 6]) # Different number of oputput prices and outputs
    @test_throws ErrorException deaprofitability([1; 2], [4 ; 5], [1; 1], [4 ; 5], Xref = [1; 2], Wref = [1; 2; 3]) # Different size in price reference sets
    @test_throws ErrorException deaprofitability([1; 2], [4 ; 5], [1; 1], [4 ; 5], Yref = [1; 2], Pref = [1; 2; 3]) # Different size in prices reference set

end
