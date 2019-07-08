# Tests for Revenue DEA Models
@testset "RevenueDEAModel" begin

    ## Test Revenue DEA Model with Cooper et al. (2007)
    # Test agains book results
    X = [3 2; 1 3; 4 6]
    Y = [3; 5; 6]
    P = [6; 6; 6]

    dearevenuecooper = dearevenue(X, Y, P, rts = :CRS)
    @test efficiency(dearevenuecooper, :Economic)   ≈ [0.9; 1; 0.6] atol = 1e-3
    @test efficiency(dearevenuecooper, :Technical)  ≈ [0.9; 1; 0.6  ] atol = 1e-3
    @test efficiency(dearevenuecooper, :Allocative) ≈ [1; 1; 1] atol = 1e-3


    ## Test Revenue DEA Model with Zofío and Prieto (2006) data.
    # Test agains results with R
    X = [5 3; 2 4; 4 2; 4 8; 7 9]
    Y = [7 4; 10 8; 8 10; 5 4; 3 6]
    P = [3 2; 3 2; 3 2; 3 2; 3 2.0]

    # Revnue CRS
    dearevenuecrs = dearevenue(X, Y, P, rts = :CRS)
    @test efficiency(dearevenuecrs) ≈ [0.4915254237;
                                    1.000;
                                    1.000;
                                    0.250;
                                    0.1735537190]
    @test efficiency(dearevenuecrs, :Technical) ≈ [0.6364;
                                1.000;
                                1.000;
                                0.250;
                                0.2609]  atol = 1e-3
    @test efficiency(dearevenuecrs, :Allocative) ≈ [0.7723970944;
                                1.000;
                                1.000;
                                1.000;
                                0.6652892562]

    # Revenue VRS
    dearevenuevrs = dearevenue(X, Y, P, rts = :VRS)
    @test efficiency(dearevenuevrs) ≈ [0.6444444444 ;
                                    1.000;
                                    1.000;
                                    0.500;
                                    0.4565217391]
    @test efficiency(dearevenuevrs, :Technical) ≈ [0.7777777778 ;
                                1.000;
                                1.000;
                                0.500 ;
                                0.6000000000]
    @test efficiency(dearevenuevrs, :Allocative) ≈ [0.8285714286;
                                1.000;
                                1.000;
                                1.000;
                                0.7608695652]

    # Check defaults
    @test efficiency(dearevenue(X, Y, P)) == efficiency(dearevenuevrs)
    @test efficiency(dearevenuevrs, :Economic) == efficiency(dearevenuevrs)

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    dearevenue_ref_eff = zeros(size(X, 1))
    dearevenue_ref_tech = zeros(size(X, 1))
    dearevenue_ref_alloc = zeros(size(X, 1))
    Xref = X[:,:]
    Yref = Y[:,:]
    Pref = P[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]
        Peval = P[i:i,:]
        Peval = Peval[:,:]

        dearevenue_ref_eff[i] = efficiency(dearevenue(Xeval, Yeval, Peval, Xref = Xref, Yref = Yref, Pref = Pref))[1]
        dearevenue_ref_tech[i] = efficiency(dearevenue(Xeval, Yeval, Peval,  Xref = Xref, Yref = Yref, Pref = Pref), :Technical)[1]
        dearevenue_ref_alloc[i] = efficiency(dearevenue(Xeval, Yeval, Peval, Xref = Xref, Yref = Yref, Pref = Pref), :Allocative)[1]

    end

    @test dearevenue_ref_eff ≈ efficiency(dearevenuevrs)
    @test dearevenue_ref_tech ≈ efficiency(dearevenuevrs, :Technical)
    @test dearevenue_ref_alloc ≈ efficiency(dearevenuevrs, :Allocative)

end
