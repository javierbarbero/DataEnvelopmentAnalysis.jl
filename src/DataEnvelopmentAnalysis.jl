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
    CostDEAModel, RevenueDEAModel, ProfitDEAModel, ProfitabilityDEAModel,
    AbstractProductivityDEAModel,
    MalmquistDEAModel,
    # Technical models
    dea, deaadd, deaddf, deagdf,
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
        # Solve nonlinear problem to display Ipopt initial message
        X = [5 3; 2 4; 4 2; 4 8; 7 9];
        Y = [7 4; 10 8; 8 10; 5 4; 3 6];
        deagdf(X, Y, 0.5, rts = :VRS)
        nothing
    end

end # module
