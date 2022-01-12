# This file contains types and structures for technical DEA model
"""
    AbstractTechnicaDEAlModel
An abstract type representing a technical DEA model.
"""
abstract type AbstractTechnicalDEAModel  <: AbstractDEAModel end

"""
    efficiency(model::AbstractTechnicalDEAModel)
Return efficiency scores of a technical DEA model.
"""
efficiency(model::AbstractTechnicalDEAModel) = model.eff

"""
    slacks(model::AbstractTechnicalDEAModel, slack::Symbol)
Return slacks of a technical DEA model.
"""
function slacks(model::AbstractTechnicalDEAModel, slack::Symbol)::Matrix

    if slack == :X
        return model.slackX
    elseif slack == :Y
        return model.slackY
    end

    throw(ArgumentError("`slack` must be :X or :Y"));

end

"""
    targets(model::AbstractTechnicalDEAModel, target::Symbol)
Return targets of a technical DEA model.
"""
function targets(model::AbstractTechnicalDEAModel, target::Symbol)::Matrix

    if target == :X
        return model.Xtarget
    elseif target == :Y
        return model.Ytarget
    end

    throw(ArgumentError("`target` must be :X or :Y"));

end

"""
    multipliers(model::AbstractTechnicalDEAModel, multiplier::Symbol)
Return multipliers (shadow prices) of a technical DEA model.
"""
function multipliers(model::AbstractTechnicalDEAModel, multiplier::Symbol)::Matrix
    if multiplier == :X
        return model.v
    elseif multiplier == :Y
        return model.u
    else
        throw(ArgumentError("`multiplier` must be :X or :Y"));
    end
end

"""
    rts(model::AbstractTechnicalDEAModel)
Return the value measuring the returns to scale of a multiplier DEA model.
"""
function rts(model::AbstractTechnicalDEAModel)::Vector
    if isdefined(model, :omega)
        return model.omega
    else
        throw(ArgumentError("`rts` only for DEA models in multiplier form."));
    end
end
