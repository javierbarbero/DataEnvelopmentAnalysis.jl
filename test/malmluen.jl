# Tests for Mamlmquist DEA Model
@testset "MalmquistLuenbergerDEAModel" begin

    ## Test Mamlmquist Luenberger DEA Model with 1 input, 1 good output and 1 bad output
    X = Array{Float64,3}(undef, 5, 1, 2)
    Y = Array{Float64,3}(undef, 5, 1, 2)
    B = Array{Float64,3}(undef, 5, 1, 2)
    
    X[:,:,1] = ones(5, 1);
    X[:,:,2] = ones(5, 1);
    
    Y[:,:,1] = [7; 5; 1; 3; 4];
    Y[:,:,2] = [8; 5.5; 2; 2; 4];
    
    B[:,:,1] = [2; 5; 3; 3; 2];
    B[:,:,2] = [1; 3; 2; 4; 1];

    # Default Malmquist Productivity Index
    mlprod = malmluen(X, Y, B)

    @test typeof(mlprod) == MalmquistLuenbergerDEAModel

    @test nobs(mlprod) == 5
    @test ninputs(mlprod) == 1
    @test noutputs(mlprod) == 1
    @test nperiods(mlprod) == 2
    @test prodchange(mlprod) ≈ [ 1.3702375781; 1.1; 1.1259778360; 0.91624569458; 1.2792042981]
    @test prodchange(mlprod, :Prod) == prodchange(mlprod)
    @test prodchange(mlprod, :EC) ≈ [1.0; 0.9625; 1.027173913; 119/144; 21/22];
    @test prodchange(mlprod, :TC) ≈ [1.3702375781;1.1428571429; 1.0961900625; 1.1087342859; 1.3401187885];

    # Test geomean is the geometric mean of TC
    mlprodbase = malmluen(X, Y, B, refperiod = :Base)
    mlprodcomparison = malmluen(X, Y, B, refperiod = :Comparison)

    @test prodchange(mlprod, :TC) == sqrt.( prodchange(mlprodbase, :TC) .* prodchange(mlprodcomparison, :TC) )

    # Test default directions
    @test prodchange(malmluen(X, Y, B, Gx = :Zeros, Gy = :Observed, Gb = :Observed)) ==  prodchange(mlprod);

    # Test other directions
    prodchange(malmluen(X, Y, B, Gx = :Zeros, Gy = :Zeros, Gb = :Zeros)) ==  [1.0; 1.0; 1.0; 1.0; 1.0];
    prodchange(malmluen(X, Y, B, Gx = :Ones, Gy = :Ones, Gb = :Ones)) ≈  [Inf; 1.0480449271555035; 1.0717826032913338; 0.9281909617845144; 1.243734296383275];
    prodchange(malmluen(X, Y, B, Gx = :Observed, Gy = :Observed, Gb = :Observed)) ≈  [1.3237752650585946; 1.0400628679223047; 1.118033988749895; 0.9045340337332908; 1.1677484162422844];
    prodchange(malmluen(X, Y, B, Gx = :Mean, Gy = :Mean, Gb = :Mean)) ≈  [1.3165611772087666; 1.0410852057520315; 1.0774756611163254; 0.9487582273837339; 1.1915664585863126];
    prodchange(malmluen(X, Y, B, Gx = ones(size(X)), Gy = ones(size(Y)), Gb = ones(size(B)))) ≈  [Inf; 1.0480449271555035; 1.0717826032913338; 0.9281909617845144; 1.243734296383275];

    ## Test Mamlmquist DEA Model with 1 input and 1 output; and 3 years
    X = Array{Float64,3}(undef, 5, 1, 3)
    Y = Array{Float64,3}(undef, 5, 1, 3)
    B = Array{Float64,3}(undef, 5, 1, 3)
    
    X[:,:,1] = ones(5, 1);
    X[:,:,2] = ones(5, 1);
    X[:,:,3] = ones(5, 1);
    
    Y[:,:,1] = [7; 5; 1; 3; 4];
    Y[:,:,2] = [8; 5.5; 2; 2; 4];
    Y[:,:,3] = [8; 5.5; 2; 2; 4];
    
    B[:,:,1] = [2; 5; 3; 3; 2];
    B[:,:,2] = [1; 3; 2; 4; 1];
    B[:,:,3] = [0.5; 2; 0.5; 3; 0.5];

    # Default Malmquist Productivity Index
    mlprod3 = malmluen(X, Y, B)

    @test nobs(mlprod3) == 5
    @test ninputs(mlprod3) == 1
    @test noutputs(mlprod3) == 1
    @test nperiods(mlprod3) == 3
    @test prodchange(mlprod3) ≈ [ 1.37024   1.22474
                                1.1       1.0
                                1.12598   1.25245
                                0.916246  1.01484
                                1.2792    1.26491] atol = 1e-5
    @test prodchange(mlprod3, :Prod) == prodchange(mlprod3)
    @test prodchange(mlprod3, :EC) ≈ [ 1.0       1.0
                                      0.9625    1.0
                                      1.02717   1.11111
                                      0.826389  0.980392
                                      0.954545  1.0] atol = 1e-5
    @test prodchange(mlprod3, :TC) ≈ [ 1.37024  1.22474
                                     1.14286  1.0
                                     1.09619  1.1272
                                     1.10873  1.03514
                                     1.34012  1.264910] atol = 1e-5

    # Print
    show(IOBuffer(), mlprod)
    show(IOBuffer(), mlprod3)

    # Test errors
    @test_throws ArgumentError prodchange(mlprod, :Error)
    @test_throws ArgumentError malmluen(X, Y, B, Gx = :Error) # Invalid inuts direction
    @test_throws ArgumentError malmluen(X, Y, B, Gy = :Error) # Invalid good outputs direction
    @test_throws ArgumentError malmluen(X, Y, B, Gb = :Error) # Invalid bad outputs direction

end
