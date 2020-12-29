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

    import StatsBase: nobs, mean


    export
    # Types
    AbstractDEAModel,
    DEAModel,
    AbstractDEAPeers, AbstractDEAPeersDMU,
    DEAPeers, DEAPeersDMU,
    AbstractTechnicalDEAModel, AbstractRadialDEAModel,
    TechnicalDEAModel,
    RadialDEAModel, AdditiveDEAModel, DirectionalDEAModel, GeneralizedDFDEAModel,
    AbstractEconomicDEAModel,
    AbstractCostDEAModel, AbstractRevenueDEAModel, AbstractProfitDEAModel, AbstractProfitabilityDEAModel,
    CostDEAModel, RevenueDEAModel, ProfitDEAModel, ProfitabilityDEAModel,
    AbstractProductivityDEAModel,
    MalmquistDEAModel,
    # All models
    nobs, ninputs, noutputs,
    # Peers
    peers, ispeer,
    # Technical models
    dea, deaadd, deaaddweights, deaddf, deagdf,
    efficiency, slacks,
    # Economic models
    deacost, dearevenue, deaprofit, deaprofitability,
    # Common technical and economic models
    targets,
    # Productivity models
    malmquist,
    nperiods, prodchange

    # Include code of functions
    include("model.jl")
    include("peers.jl")

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
