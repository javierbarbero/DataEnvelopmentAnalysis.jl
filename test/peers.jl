# Tests for DEAPeers and DEAPeersDMU
@testset "DEAPeers" begin

    # Test agains results in R
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    P = peers(dea(X, Y, orient = :Input, rts = :CRS))

    @test typeof(P) == DEAPeers
    @test eltype(P) == DEAPeersDMU

    @test size(P) == (11,)
    @test ispeer(P, 1, 1) == true
    @test ispeer(P, 1, 2) == false
    @test ispeer(P, 2, 1) == false
    @test ispeer(P, 2, 3) == false

    @test convert(Matrix, P) ≈
    [1.000000000  0  0 0.0000000000  0  0 0.00000000000  0  0   0   0
    0.000000000  0  0 0.4249783174  0  0 0.10928013877  0  0   0   0
    1.134321653  0  0 0.4380053908  0  0 0.00000000000  0  0   0   0
    0.000000000  0  0 1.0000000000  0  0 0.00000000000  0  0   0   0
    0.000000000  0  0 0.2573807721  0  0 0.04844814534  0  0   0   0
    0.000000000  0  0 0.0000000000  0  0 0.33333333333  0  0   0   0
    0.000000000  0  0 0.0000000000  0  0 1.00000000000  0  0   0   0
    0.000000000  0  0 1.0348650979  0  0 0.11457435013  0  0   0   0
    0.000000000  0  0 0.0000000000  0  0 1.14814814815  0  0   0   0
    0.000000000  0  0 0.4905660377  0  0 0.49056603774  0  0   0   0
    0.000000000  0  0 0.0000000000  0  0 0.00000000000  0  0   0   1]

    # Test using DMU names
    firms = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; "I"; "J"; "L"]

    Pn = peers(dea(X, Y, orient = :Input, rts = :CRS, names = firms))

    @test typeof(P) == DEAPeers
    @test eltype(P) == DEAPeersDMU

    @test size(P) == (11,)
    @test length(P) == 11
    @test firstindex(P) == 1
    @test lastindex(P) == 11

    @test ispeer(Pn, "A", "A") == true
    @test ispeer(Pn, "A", "B") == false
    @test ispeer(Pn, "B", "A") == false
    @test ispeer(Pn, "B", "C") == false

    # Test peers of a specific DMU
    P2 = P[2]

    @test typeof(P2) == DEAPeersDMU
    @test eltype(P2) == Tuple{Tuple{Int64,String},Float64}

    @test size(P2) == (2,)
    @test length(P2) == 2
    @test firstindex(P2) == 1
    @test lastindex(P2) == 2

    @test ispeer(P2, 1) == false
    @test ispeer(P2, 2) == false
    @test ispeer(P2, 3) == false
    @test ispeer(P2, 4) == true
    @test ispeer(P2, 7) == true

    # Test peers of a specific DMU using names
    P2n = Pn[2]

    @test ispeer(P2n, "A") == false
    @test ispeer(P2n, "B") == false
    @test ispeer(P2n, "C") == false
    @test ispeer(P2n, "D") == true
    @test ispeer(P2n, "G") == true

    # Test elements of DEAPeersDMU
    @test P2n[1] == ((4, "D"), 0.424978317432784)
    @test P2n[2] == ((7, "G"), 0.10928013876843023)

    # Sum function
    @test sum(P, dims = 2) ≈ [ 1.0
            0.5342584562012143
            1.5723270440251573
            1.0
            0.30582891748675245
            0.3333333333333333
            1.0
            1.1494394480719479
            1.1481481481481481
            0.9811320754716981
            1.0]

    @test sum(peers(dea(X, Y, rts = :VRS)), dims = 2) ≈ ones(11, 1)

    # Test references names
    Pnamesref = peers(dea(X, Y), namesref =["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"])
    @test Pnamesref.dmunamesref == ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]
    @test Pnamesref[11][1][1][2] == "K"

    # Test conversions
    @test typeof(convert(Matrix, P)) == Array{Float64,2}
    @test typeof(convert(SparseMatrixCSC, P)) == SparseMatrixCSC{Float64,Int64}

    # Print
    show(IOBuffer(), P)
    show(IOBuffer(), P2)

    # Test errors
    @test_throws ErrorException peers(dea(X, Y, rts = :VRS), namesref = ["A"]) #  Length of references names different to number of references DMUs
    @test_throws ErrorException ispeer(Pn, "W", "A") # First name does not exists
    @test_throws ErrorException ispeer(Pn, "A", "W") # Second name does not exists

    Prepeated = peers(dea([1; 1; 1; 1; 1], [1 ; 2 ;3; 4; 5], orient = :Output, names = ["A", "B", "B", "C", "C"]))
    @test_throws ErrorException ispeer(Prepeated, "C", "A") # First name is repeated
    @test_throws ErrorException ispeer(Prepeated, "A", "C") # Second name does not exists

    @test_throws BoundsError ispeer(P2, 0) # BoundsErrror

    # Test struct defaults and errors
    struct wrongDEAmodel <: AbstractDEAModel
        n::Int64
        eff::Vector
    end

    @test names(wrongDEAmodel(3, [1; 2; 3])) == ["1"; "2"; "3"] # Default names if dmunames not in struct
    @test_throws ErrorException peers(wrongDEAmodel(3, [1; 2; 3])) # Model does not have info on peers

end
