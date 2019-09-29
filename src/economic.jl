# This file contains types and structures for economic DEA model
"""
    AbstractEconomicDEAlModel
An abstract type representing an economic DEA model.
"""
abstract type AbstractEconomicDEAModel <: AbstractTechnicalDEAModel end

"""
    efficiency(model::AbstractEconomicDEAModel)
Return efficiency scores of an economic DEA model.

# Optional Arguments
- `type=Economic`: type of efficiency scores to return.

Type specification:
- `:Economic`: returns economic efficiency of the model.
- `:CRS`: returns technical efficiency under constant returns to scale.
- `:VRS`: returns technical efficiency under variable returns to scale.
- `:Scale`: returns scale efficiency.
- `:Allocative`: returns allocative efficiency.

# Examples
```jldoctest
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9.0];
julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];
julia> W = [2 1; 2 1; 2 1; 2 1; 2 1.0];
julia> P = [3 2; 3 2; 3 2; 3 2; 3 2.0];
julia> profitbl = deaprofitability(X, Y, W, P)
julia> efficiency(profitbl)
5-element Array{Float64,1}:
 0.38795983677810825
 0.9999999082180037
 0.7652173234322985
 0.249999984621284
 0.15879016408241559

julia> efficiency(profitbl, :Allocative)
5-element Array{Float64,1}:
 0.6096511887087087
 0.9999999387832732
 0.765217344390068
 0.9999999256504388
 0.6086956149040872
```
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

    error(typeof(model), " has no efficiency type ", string(type))


end
