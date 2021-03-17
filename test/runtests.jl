using DataEnvelopmentAnalysis
using LinearAlgebra
using SparseArrays
using Test

@testset "DataEnvelopmentAnalysis" begin

    include("dea.jl")
    include("deaadd.jl")
    include("deaddf.jl")
    include("deagdf.jl")
    include("dearussell.jl")
    include("deaerg.jl")
    include("deamddf.jl")

    include("econoptim.jl")
    include("deacost.jl")
    include("dearevenue.jl")
    include("deaprofit.jl")
    include("deaprofitability.jl")

    include("malmquist.jl")

    include("peers.jl")

end
