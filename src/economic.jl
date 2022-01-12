# This file contains types and structures for economic DEA models
"""
    AbstractEconomicDEAlModel
An abstract type representing an economic DEA model.
"""
abstract type AbstractEconomicDEAModel <: AbstractDEAModel end

"""
    AbstractCostDEAModel
An abstract type representing a cost DEA model.
"""
abstract type AbstractCostDEAModel <: AbstractEconomicDEAModel end

"""
    AbstractRevenueDEAModel
An abstract type representing a revenue DEA model.
"""
abstract type AbstractRevenueDEAModel <: AbstractEconomicDEAModel end

"""
    AbstractProfitDEAModel
An abstract type representing a revenue DEA model.
"""
abstract type AbstractProfitDEAModel <: AbstractEconomicDEAModel end

"""
    AbstractProfitabilityDEAModel
An abstract type representing a revenue DEA model.
"""
abstract type AbstractProfitabilityDEAModel <: AbstractEconomicDEAModel end

"""
    efficiency(model::AbstractEconomicDEAModel)
Return efficiency scores of an economic DEA model.

# Optional Arguments
- `type=Economic`: type of efficiency scores to return.

Type specification:
- `:Economic`: returns economic efficiency of the model.
- `:Technical`: returns technical efficiency.
- `:Allocative`: returns allocative efficiency.

Some models also allow these types:
- `:CRS`: returns technical efficiency under constant returns to scale.
- `:VRS`: returns technical efficiency under variable returns to scale.
- `:Scale`: returns scale efficiency.
"""
function efficiency(model::AbstractEconomicDEAModel, type::Symbol = :Economic)::Vector

    if type == :Economic
        return model.eff
    end

    if type == :Technical
        if isdefined(model, :techeff)
            return model.techeff
        end
    end

    if type == :CRS
        if isdefined(model, :crseff)
            return model.crseff
        end
    end

    if type == :VRS
        if isdefined(model, :vrseff)
            return model.vrseff
        end
    end

    if type == :Scale
        if isdefined(model, :scaleff)
            return model.scaleff
        end
    end

    if type == :Allocative
        if isdefined(model, :alloceff)
            return model.alloceff
        end
    end

    throw(ArgumentError("$(typeof(model)) has no efficiency $(type)"));

end

"""
    targets(model::AbstractEconomicDEAModel, target::Symbol)
Return targets of an economic DEA model.
"""
function targets(model::AbstractEconomicDEAModel, target::Symbol)::Matrix

    if target == :X
        return model.Xtarget
    elseif target == :Y
        return model.Ytarget
    end

    throw(ArgumentError("`target` must be :X or :Y"));

end

"""
    normfactor(model::AbstractEconomicDEAModel)
Return the normalization factor of an economic DEA model.
"""
function normfactor(model::AbstractEconomicDEAModel)::Vector
    if isdefined(model, :normalization)
        return model.normalization
    else
        throw(ArgumentError("$(typeof(model)) has no normalization factor"));
    end
end

"""
    ismonetary(model::AbstractEconomicDEAModel)
Indicate whether inefficiency is in monetary units.
"""
function ismonetary(model::AbstractEconomicDEAModel)::Bool
    if isdefined(model, :monetary)
        return model.monetary
    else
        throw(ArgumentError("$(typeof(model)) has no monetary identifier"));
    end
end
