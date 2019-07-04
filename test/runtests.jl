using DataEnvelopmentAnalysis
using LinearAlgebra
using Test

@testset "DataEnvelopmentAnalysis" begin

    include("dea.jl")
    include("deaadd.jl")
    include("deagdf.jl")

end
