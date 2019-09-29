# Tests for Profit DEA Models
@testset "ProfitDEAModel" begin

    ## Test Profit DEA Model
    # Test with Zofio, Pastor and Aparicio (2013) data
    X = [1 1; 1 1; 0.75 1.5; 0.5 2; 0.5 2; 2 2; 2.75 3.5; 1.375 1.75]
    Y = [1 11; 5 3; 5 5; 2 9; 4 5; 4 2; 3 3; 4.5 3.5]
    P = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1]
    W = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1]

    GxGydollar = 1 ./ (sum(P, dims = 2) + sum(W, dims = 2))
    GxGydollar = repeat(GxGydollar, 1, 2)

    deaprofitdollar = deaprofit(X, Y, W, P, GxGydollar, GxGydollar)
    @test efficiency(deaprofitdollar, :Economic)   ≈ [2; 2; 0; 2; 2; 8; 12; 4] atol = 1e-3
    @test efficiency(deaprofitdollar, :Technical)  ≈ [0; 0; 0; 0; 0; 6; 12; 3] atol = 1e-3
    @test efficiency(deaprofitdollar, :Allocative) ≈ [2; 2; 0; 2; 2; 2; 0; 1] atol = 1e-3

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deaprofit_ref_eff = zeros(size(X, 1))
    deaprofit_ref_tech = zeros(size(X, 1))
    deaprofit_ref_alloc = zeros(size(X, 1))
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
        Gxeval = GxGydollar[i:i,:]
        Gyeval = GxGydollar[i:i,:]

        deaprofit_ref_eff[i] = efficiency(deaprofit(Xeval, Yeval, Weval, Peval, Gxeval, Gyeval, Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref))[1]
        deaprofit_ref_tech[i] = efficiency(deaprofit(Xeval, Yeval, Weval, Peval,  Gxeval, Gyeval, Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref), :Technical)[1]
        deaprofit_ref_alloc[i] = efficiency(deaprofit(Xeval, Yeval, Weval, Peval, Gxeval, Gyeval, Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref), :Allocative)[1]

    end

    @test deaprofit_ref_eff ≈ efficiency(deaprofitdollar)
    @test deaprofit_ref_tech ≈ efficiency(deaprofitdollar, :Technical)
    @test deaprofit_ref_alloc ≈ efficiency(deaprofitdollar, :Allocative)

    # Print
    show(IOBuffer(), deaprofitdollar)

end
