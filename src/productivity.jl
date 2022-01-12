# This file contains types and structures for productivity DEA model
"""
    AbstractProductivityDEAlModel
An abstract type representing a productivity DEA model.
"""
abstract type AbstractProductivityDEAModel <: AbstractDEAModel end

"""
    nperiods(model::AbstractProductivityDEAModel)
Return number of time periods of a productivity DEA model.
"""
nperiods(model::AbstractProductivityDEAModel) = model.periods

"""
    prodchange(model::AbstractProductivityDEAModel)
Return productivity change of a productivity change DEA model.
"""
function prodchange(model::AbstractProductivityDEAModel, type::Symbol = :Prod)::Matrix

    if type == :Prod
        return model.Prod
    end

    if type == :EC
        if isdefined(model, :EC)
            return model.EC
        end
    end

    if type == :TC
        if isdefined(model, :TC)
            return model.TC
        end
    end

    throw(ArgumentError("$(typeof(model)) has no productivity change component $(type)"));

end
