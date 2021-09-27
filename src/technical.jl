# This file contains types and structures for technical DEA model
"""
    AbstractTechnicaDEAlModel
An abstract type representing a technical DEA model.
"""
abstract type AbstractTechnicalDEAModel  <: AbstractDEAModel end

"""
    efficiency(model::AbstractTechnicalDEAModel)
Return efficiency scores of a technical DEA model.

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> efficiency(deaio)
11-element Vector{Float64}:
 1.0
 0.6222896790980051
 0.8198562443845464
 1.0
 ⋮
 0.7576690895651103
 0.8201058201058201
 0.49056603773584917
 1.0
```
"""
efficiency(model::AbstractTechnicalDEAModel) = model.eff

"""
    slacks(model::AbstractTechnicalDEAModel, slack::Symbol)
Return slacks of a technical DEA model.
# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> slacks(deaio, :X)
11×2 Matrix{Float64}:
  0.0          0.0
 -4.41868e-15  0.0
  0.0          8.17926e-15
 -8.03397e-16  0.0
  ⋮            
  1.60679e-15  0.0
  1.64021      0.0
  9.68683e-15  0.0
  0.0          4.0

julia> slacks(deaio, :Y)
11×1 Matrix{Float64}:
 0.0
 0.0
 0.0
 0.0
 ⋮
 0.0
 0.0
 0.0
 0.0
```
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
# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> targets(deaio, :X)
11×2 Matrix{Float64}:
  5.0      13.0
  9.95663   7.46748
 13.1177   21.3163
 17.0      15.0
  ⋮        
 20.4571   16.6687
 28.7037   11.4815
 20.6038   12.2642
  5.0      13.0

julia> targets(deaio, :Y)
11×1 Matrix{Float64}:
 12.0
 14.0
 25.0
 26.0
  ⋮
 30.0
 31.0
 26.0
 12.0
```
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
