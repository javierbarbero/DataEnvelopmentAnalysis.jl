# Tests for Cost DEA Models
@testset "CostDEAModel" begin

    ## Test Cost DEA Model with Cooper et al. (2007)
    # Test agains book results
    X = [3 2; 1 3; 4 6]
    Y = [3; 5; 6]
    W = [4 2; 4 2; 4 2]

    deacostcooper = deacost(X, Y, W, rts = :CRS)
    @test efficiency(deacostcooper, :Economic)   ≈ [0.375; 1; 0.429] atol = 1e-3
    @test efficiency(deacostcooper, :Technical)  ≈ [0.9  ; 1; 0.6  ] atol = 1e-3
    @test efficiency(deacostcooper, :Allocative) ≈ [0.417; 1; 0.714] atol = 1e-3


    ## Test Cost DEA Model with Zofío and Prieto (2006) data.
    # Test agains results with R
    X = [5 3; 2 4; 4 2; 4 8; 7 9]
    Y = [7 4; 10 8; 8 10; 5 4; 3 6]
    W = [2 1; 2 1; 2 1; 2 1; 2 1.0]

    # Cost CRS
    deacostcrs = deacost(X, Y, W, rts = :CRS)
    @test efficiency(deacostcrs) ≈ [0.4307692308;
                                    1.000;
                                    1.000;
                                    0.250;
                                    0.2608695652]
    @test efficiency(deacostcrs, :Technical) ≈ [0.6364;
                                1.000;
                                1.000;
                                0.250;
                                0.2609]  atol = 1e-3
    @test efficiency(deacostcrs, :Allocative) ≈ [0.6769230769;
                                1.000;
                                1.000;
                                1.000;
                                1.000]

    # Cost VRS
    deacostvrs = deacost(X, Y, W, rts = :VRS)
    @test efficiency(deacostvrs) ≈ [0.6153846154;
                                    1.000;
                                    1.000;
                                    0.500;
                                    0.3478260870]
    @test efficiency(deacostvrs, :Technical) ≈ [0.750;
                                1.000;
                                1.000;
                                0.500 ;
                                0.375]  atol = 1e-3
    @test efficiency(deacostvrs, :Allocative) ≈ [0.8205128205;
                                1.000;
                                1.000;
                                1.000;
                                0.9275362319]

    # Check defaults
    @test efficiency(deacost(X, Y, W)) == efficiency(deacostvrs)
    @test efficiency(deacostvrs, :Economic) == efficiency(deacostvrs)

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deacost_ref_eff = zeros(size(X, 1))
    deacost_ref_tech = zeros(size(X, 1))
    deacost_ref_alloc = zeros(size(X, 1))
    Xref = X[:,:]
    Yref = Y[:,:]
    Wref = W[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]
        Weval = W[i:i,:]
        Weval = Weval[:,:]

        deacost_ref_eff[i] = efficiency(deacost(Xeval, Yeval, Weval, Xref = Xref, Yref = Yref, Wref = Wref))[1]
        deacost_ref_tech[i] = efficiency(deacost(Xeval, Yeval, Weval,  Xref = Xref, Yref = Yref, Wref = Wref), :Technical)[1]
        deacost_ref_alloc[i] = efficiency(deacost(Xeval, Yeval, Weval, Xref = Xref, Yref = Yref, Wref = Wref), :Allocative)[1]

    end

    @test deacost_ref_eff ≈ efficiency(deacostvrs)
    @test deacost_ref_tech ≈ efficiency(deacostvrs, :Technical)
    @test deacost_ref_alloc ≈ efficiency(deacostvrs, :Allocative)

    # Print
    show(IOBuffer(), deacostcooper)

    # Test errors
    @test_throws ErrorException deacost([1; 2 ; 3], [4 ; 5], [1; 1; 1]) #  Different number of observations
    @test_throws ErrorException deacost([1; 2], [4 ; 5], [1; 1], Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws ErrorException deacost([1 1; 2 2], [4 4; 5 5], [1 1; 2 2], Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws ErrorException deacost([1 1; 2 2], [4 4; 5 5], [4 4; 5 5], Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ErrorException deacost([1; 2; 3], [4; 5; 6], [1; 2; 3], rts = :Error) # Invalid returns to scale
    @test_throws ErrorException deacost([1; 2; 3], [4; 5; 6], [1; 2; 3; 4]) # Different number of observation in prices
    @test_throws ErrorException deacost([1 1; 2 2; 3 3 ], [4; 5; 6], [1 1 1; 2 2 2; 3 3 3]) # Different number of input prices and inputs
    @test_throws ErrorException deacost([1; 2], [4 ; 5], [1; 1], Xref = [1; 2], Wref = [1; 2; 3]) # Different size in price reference sets

end
