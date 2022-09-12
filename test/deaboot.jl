# Tests for BootstrapRadial DEA Models
@testset "BootstrapRadialDEAModel" begin

    ## Test Bootstrap Radial DEA Model
    X = [2, 4, 3, 5, 6]
    Y = [1, 2, 3, 4, 5]

    ioboot = deaboot(X, Y, orient = :Input, rts = :VRS, rng = StableRNG(1234567))

    @test efficiency(ioboot) ≈ [0.936384; 0.606179; 0.940476; 0.869816; 0.907370] atol = 1e-5
    @test ioboot.effref ≈ [1.0; 0.625; 1.0; 0.9; 1.0] atol = 1e-5
    @test bandwidth(ioboot) ≈ 0.1341004 atol = 1e-5
    @test confint(ioboot) ≈ [0.75  1.0; 0.535644  0.625; 0.816986  1.0; 0.706839  0.9; 0.734459  1.0] atol = 1e-5

    ooboot = deaboot(X, Y, orient = :Output, rts = :VRS, rng = StableRNG(1234567))

    @test efficiency(ooboot) ≈ [1.096844; 1.874885; 1.075607; 1.108465; 1.05909] atol = 1e-5
    @test ooboot.effref ≈ [1; 11/6; 1; 13/12; 1.0] atol = 1e-5
    @test bandwidth(ooboot) ≈ 0.1398634 atol = 1e-5
    @test confint(ooboot) ≈ [1.0  1.44558; 1.83333  2.17557; 1.0  1.32604; 1.08333  1.26799; 1.0  1.28106] atol = 1e-5

end
