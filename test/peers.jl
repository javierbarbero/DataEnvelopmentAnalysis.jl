# Tests for DEAPeers and DEAPeersDMU
@testset "DEAPeers" begin

    # Test agains results in R
    X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17]
    Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12]

    P = peers(dea(X, Y, orient = :Input, rts = :CRS))

    @test typeof(P) == DEAPeers
    @test eltype(P) == DEAPeersDMU
    @test IndexStyle(P) == IndexLinear()

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

    @test ispeer(Pn, "A", "A") == true
    @test ispeer(Pn, "A", "B") == false
    @test ispeer(Pn, "B", "A") == false
    @test ispeer(Pn, "B", "C") == false

    # Test peers of a specific DMU
    P2 = P[2]

    @test typeof(P2) == DEAPeersDMU
    @test eltype(P2) == Tuple{Tuple{Int64,String},Float64}
    @test IndexStyle(P) == IndexLinear()

    @test size(P2) == (2,)
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

    # Test conversions
    @test typeof(convert(Matrix, P)) == Array{Float64,2}
    #@test typeof(convert(SparseMatrixCSC, P)) == SparseMatrixCSC{Float64,Int64}

    # Print
    show(IOBuffer(), P)
    show(IOBuffer(), P2)

end
