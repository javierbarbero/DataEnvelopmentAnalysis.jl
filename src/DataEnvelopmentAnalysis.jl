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
    TechnicalDEAModel,
    RadialDEAModel, AdditiveDEAModel, DirectionalDEAModel, GeneralizedDFDEAModel,
    AbstractEconomicDEAModel,
    AbstractCostDEAModel, AbstractRevenueDEAModel, AbstractProfitDEAModel, AbstractProfitabilityDEAModel,
    CostDEAModel, RevenueDEAModel, ProfitDEAModel, ProfitabilityDEAModel,
    AbstractProductivityDEAModel,
    MalmquistDEAModel,
    # Technical models
    dea, deaadd, deaaddweights, deaddf, deagdf,
    efficiency, slacks,
    nobs, ninputs, noutputs, peers,
    # Economic models
    deacost, dearevenue, deaprofit, deaprofitability,
    # Productivity models
    malmquist,
    nperiods, prodchange

    # Include code of functions
    include("technical.jl")
    include("dea.jl")
    include("deaadd.jl")
    include("deaddf.jl")
    include("deagdf.jl")

    include("economic.jl")
    include("deacost.jl")
    include("dearevenue.jl")
    include("deaprofit.jl")
    include("deaprofitability.jl")

    include("productivity.jl")
    include("malmquist.jl")

    function __init__()

        nothing
    end

end # module
