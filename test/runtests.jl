using DataEnvelopmentAnalysis
using Distributions
using LinearAlgebra
using Random
using SparseArrays
using StableRNGs
using Test

@testset "DataEnvelopmentAnalysis" begin

    include("dea.jl")
    include("deam.jl")
    include("deaadd.jl")
    include("deaddf.jl")
    include("deaddfm.jl")
    include("deagdf.jl")
    include("dearussell.jl")
    include("deaerg.jl")
    include("deamddf.jl")
    include("deaholder.jl")
    include("dearddf.jl")

    include("econoptim.jl")
    include("deacost.jl")
    include("dearevenue.jl")
    include("deaprofit.jl")
    include("deaprofitability.jl")

    include("deaenv.jl")

    include("malmquist.jl")
    include("malmluen.jl")

    include("deaboot.jl")
    include("deartstest.jl")

    include("peers.jl")

    include("deabigdata.jl")

end
