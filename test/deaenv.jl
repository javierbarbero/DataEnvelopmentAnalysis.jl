# Tests for Environmental DEA Models
@testset "EnvironmentalDEAModel" begin

    ## Test Environmental DEA Model
    X = ones(5, 1);
    Y = [7; 5; 1; 3; 4];
    B = [2; 5; 3; 3; 2];

    # Environmental DDF model
    deaenv1 = deaenv(X, Y, B)

    @test typeof(deaenv1) == EnvironmentalDEAModel
    @test isenvironmental(deaenv1) == true;

    @test nobs(deaenv1) == 5
    @test ninputs(deaenv1) == 1
    @test noutputs(deaenv1) == 1
    @test efficiency(deaenv1) ≈ [0.0; 0.4000000000000002; 0.8260869565217392; 0.5555555555555556; 0.27272727272727276]
    @test convert(Matrix, peers(deaenv1)) ≈
                    [ 1.0       0.0  0.0  0.0  0.0;
                    1.0       0.0  0.0  0.0  0.0;
                    0.26087   0.0  0.0  0.0  0.0;
                    0.666667  0.0  0.0  0.0  0.0;
                    0.727273  0.0  0.0  0.0  0.0] atol = 1e-5
    @test slacks(deaenv1, :X) ≈ [0.0; 0.0; 0.7391304347826089; 1/3; 3/11]
    @test slacks(deaenv1, :Y) ≈ zeros(5) atol = 1e-10
    @test slacks(deaenv1, :B) ≈ [0.0; 1.0; 0.0; 0.0; 0.0]

    @test efficiency(deaenv(targets(deaenv1, :X), targets(deaenv1, :Y), targets(deaenv1, :B), slack = false)) ≈ zeros(5) atol = 1e-10
    @test efficiency(deaadd([targets(deaenv1, :X) targets(deaenv1, :B)], targets(deaenv1, :Y))) ≈ zeros(5) atol=1e-14

    @test peersmatrix(deaenv1) == deaenv1.lambda

    ## Default directions
    @test efficiency(deaenv(X, Y, B, Gx = :Zeros, Gy = :Observed, Gb = :Observed)) ≈ efficiency(deaenv1)

    ## Test other directions
    @test efficiency(deaenv(X, Y, B, Gx = :Zeros, Gy = :Zeros, Gb = :Zeros)) ≈ zeros(5)
    @test efficiency(deaenv(X, Y, B, Gx = :Ones, Gy = :Ones, Gb = :Ones)) ≈ [0; 1/4; 3/4; 0.5; 3/8]
    @test efficiency(deaenv(X, Y, B, Gx = :Observed, Gy = :Observed, Gb = :Observed)) ≈ [0; 1/6; 3/4; 2/5; 3/11] 
    @test efficiency(deaenv(X, Y, B, Gx = :Mean, Gy = :Mean, Gb = :Mean)) ≈ [0; 2/11; 6/11; 4/11; 0.20689655172413798] atol = 1e-10
    @test efficiency(deaenv(X, Y, B, Gx = X, Gy = Y, Gb = B)) ≈ [0; 1/6; 3/4; 2/5; 3/11] # Custom direction

    ## Test no slacks
    deaenv1noslacks = deaenv(X, Y, B, slack = false)
    @test efficiency(deaenv1noslacks) == efficiency(deaenv1)
    @test isempty(slacks(deaenv1noslacks, :X)) == 1
    @test isempty(slacks(deaenv1noslacks, :Y)) == 1

    @test efficiency(deaenv(targets(deaenv1noslacks, :X), targets(deaenv1noslacks, :Y), targets(deaenv1noslacks, :B), slack = false)) ≈ zeros(5) atol = 1e-10
    @test efficiency(deaadd([targets(deaenv1noslacks, :X) targets(deaenv1noslacks, :B)], targets(deaenv1noslacks, :Y))) != zeros(5) # Different as there is no slacks in first model

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deaenv_ref_eff = zeros(size(X, 1))
    deaenv_ref_eff = zeros(size(X, 1))

    deaenv_ref_slackX = zeros(size(X))
    deaenv_ref_slackY = zeros(size(Y))
    deaenv_ref_slackB = zeros(size(Y))

    Xref = X[:,:]
    Yref = Y[:,:]
    Bref = B[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]
        Beval = B[i:i,:]
        Beval = Beval[:,:]

        deaenv_ref_eff[i] = efficiency(deaenv(Xeval, Yeval, Beval, Xref = Xref, Yref = Yref, Bref = Bref))[1]

        deaenv_ref_slackX[i,:] = slacks(deaenv(Xeval, Yeval, Beval, Xref = Xref, Yref = Yref, Bref = Bref), :X)
        deaenv_ref_slackY[i,:] = slacks(deaenv(Xeval, Yeval, Beval, Xref = Xref, Yref = Yref, Bref = Bref), :Y)
        deaenv_ref_slackB[i,:] = slacks(deaenv(Xeval, Yeval, Beval, Xref = Xref, Yref = Yref, Bref = Bref), :B)
    end

    @test deaenv_ref_eff ≈ efficiency(deaenv1)

    @test deaenv_ref_slackX ≈ slacks(deaenv1, :X) atol=1e-14
    @test deaenv_ref_slackY ≈ slacks(deaenv1, :Y) atol=1e-14
    @test deaenv_ref_slackB ≈ slacks(deaenv1, :B) atol=1e-14

    # Print
    show(IOBuffer(), deaenv1)
    show(IOBuffer(), deaenv1noslacks)

    # Test errors
    @test_throws ArgumentError deaenv([1; 2; 3], [4; 5; 6], [7; 8; 9], rts = :VRS) # No VRS
    @test_throws ArgumentError deaenv([1; 2; 3], [4; 5; 6], [7; 8; 9], rts = :Error) # Invalid returns to scale

    @test_throws ArgumentError deaenv([1; 2; 3], [1; 2; 3], [7; 8; 9], Gx = :Error) # Invalid inuts direction
    @test_throws ArgumentError deaenv([1; 2; 3], [1; 2; 3], [7; 8; 9], Gy = :Error) # Invalid outputs direction
    @test_throws ArgumentError deaenv([1; 2; 3], [1; 2; 3], [7; 8; 9], Gb = :Error) # Invalid bad outputs direction

end
