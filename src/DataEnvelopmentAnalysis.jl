module DataEnvelopmentAnalysis

    """
        DataEnvelopmentAnalysis
    A Julia package for efficiency and productivity measurement using Data Envelopment Analysis (DEA).
    [DataEnvelopmentAnalysis repository](https://github.com/javierbarbero/DataEnvelopmentAnalysis.jl).
    """

    using JuMP
    using GLPK
    using Ipopt
    using SparseArrays
    using LinearAlgebra

    using Printf: @sprintf
    using Statistics: std
    using StatsBase: CoefTable

    import StatsBase: nobs


    export
    # Types
    AbstractTechnicalDEAModel, AbstractRadialDEAModel,
    TechnicalDEAModel, RadialDEAModel,
    AbstractEconomicDEAModel,
    ProfitabilityDEAModel,
    # Technical models
    dea, deaadd, deagdf,
    efficiency,
    nobs, ninputs, noutputs, peers,
    # Economic models
    deaprofitability

    # Include code of functions
    include("technical.jl")
    include("dea.jl")
    include("deaadd.jl")
    include("deagdf.jl")
    include("economic.jl")
    include("deaprofitability.jl")

end # module
