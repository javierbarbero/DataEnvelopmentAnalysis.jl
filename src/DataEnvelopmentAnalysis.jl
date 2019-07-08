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
    CostDEAModel, ProfitabilityDEAModel,
    # Technical models
    dea, deaadd, deagdf,
    efficiency,
    nobs, ninputs, noutputs, peers,
    # Economic models
    deacost, deaprofitability

    # Include code of functions
    include("technical.jl")
    include("dea.jl")
    include("deaadd.jl")
    include("deagdf.jl")
    include("economic.jl")
    include("deacost.jl")
    include("deaprofitability.jl")

    function __init__()
        # Solve nonlinear problem to display Ipopt initial message
        X = [5 3; 2 4; 4 2; 4 8; 7 9];
        Y = [7 4; 10 8; 8 10; 5 4; 3 6];
        deagdf(X, Y, 0.5, rts = :VRS)
        nothing
    end

end # module
