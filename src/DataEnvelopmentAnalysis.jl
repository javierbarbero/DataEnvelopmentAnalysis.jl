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
    using InvertedIndices 
    using ProgressMeter

    using Printf: @sprintf
    using Statistics: std, quantile!
    using StatsBase: CoefTable

    import StatsBase: nobs, mean


    export
    # optimizer
    DEAOptimizer, 
    newdeamodel,

    # Types
    AbstractDEAModel,
    AbstractDEAPeers, AbstractDEAPeersDMU,
    DEAPeers, DEAPeersDMU,
    AbstractTechnicalDEAModel, AbstractRadialDEAModel, AbstractRadialMultiplierDEAModel,
    RadialDEAModel, RadialMultiplierDEAModel, AdditiveDEAModel, 
    DirectionalDEAModel, DirectionalMultiplierDEAModel,
    GeneralizedDFDEAModel, 
    RussellDEAModel, EnhancedRussellGraphDEAModel, ModifiedDDFDEAModel, 
    AbstractHolderDEAModel, HolderL1DEAModel, HolderL2DEAModel, HolderLInfDEAModel,
    ReverseDDFDEAModel,
    AbstractEconomicDEAModel,
    AbstractCostDEAModel, AbstractRevenueDEAModel, AbstractProfitDEAModel, AbstractProfitabilityDEAModel,
    CostDEAModel, RevenueDEAModel, ProfitDEAModel, ProfitabilityDEAModel,
    AbstractProductivityDEAModel,
    MalmquistDEAModel,
    # All models
    nobs, ninputs, noutputs,
    # Peers
    peers, peersmatrix, ispeer,
    # Technical models
    dea, deam, deaadd, deaaddweights, deaddf, deaddfm, deagdf, 
    dearussell, deaerg, deamddf, deaholder, dearddf,
    deabigdata,
    efficiency, slacks, multipliers, rts,
    # Economic models
    deamincost, deamaxrevenue, deamaxprofit,
    deacost, dearevenue, deaprofit, deaprofitability, 
    normfactor, ismonetary,
    # Common technical and economic models
    targets,
    # Productivity models
    malmquist,
    nperiods, prodchange



    # Include code of functions
    include("optimizer.jl")

    include("model.jl")
    include("peers.jl")

    include("technical.jl")
    include("dea.jl")
    include("deam.jl")
    include("deaadd.jl")
    include("deaddfm.jl")
    include("deabigdata.jl")
    include("deaddf.jl")
    include("deagdf.jl")
    include("dearussell.jl")
    include("deaerg.jl")
    include("deamddf.jl")
    include("deaholder.jl")
    include("dearddf.jl")

    include("economic.jl")
    include("econoptim.jl")
    include("deacost.jl")
    include("dearevenue.jl")
    include("deaprofit.jl")
    include("deaprofitability.jl")

    include("productivity.jl")
    include("malmquist.jl")

    include("progressbarmeter.jl")
    
    function __init__()

        nothing
    end

end # module
