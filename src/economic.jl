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

# Examples
```jldoctest
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9.0];
julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];
julia> W = [2 1; 2 1; 2 1; 2 1; 2 1.0];
julia> P = [3 2; 3 2; 3 2; 3 2; 3 2.0];
julia> profitbl = deaprofitability(X, Y, W, P)
julia> efficiency(profitbl)
5-element Vector{Float64}:
 0.38795983677810825
 0.9999999082180037
 0.7652173234322985
 0.24999998462128403
 0.15879016408241559

julia> efficiency(profitbl, :Allocative)
5-element Vector{Float64}:
 0.6096511887087086
 0.9999999387832732
 0.765217344390068
 0.9999999256504392
 0.608695614904087
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

    throw(ArgumentError("$(typeof(model)) has no efficiency $(type)"));

end

"""
    targets(model::AbstractEconomicDEAModel, target::Symbol)
Return targets of an economic DEA model.
# Examples
```jldoctest
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9.0];
julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];
julia> W = [2 1; 2 1; 2 1; 2 1; 2 1.0];
julia> P = [3 2; 3 2; 3 2; 3 2; 3 2.0];
julia> profit = deaprofit(X, Y, W, P, Gx = :Monetary, Gy = :Monetary);
julia> targets(profit, :X)
5×2 Matrix{Float64}:
 2.0  4.0
 2.0  4.0
 2.0  4.0
 2.0  4.0
 2.0  4.0

julia> targets(profit, :Y)
5×2 Matrix{Float64}:
 10.0  8.0
 10.0  8.0
 10.0  8.0
 10.0  8.0
 10.0  8.0

```
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

# Examples
```jldoctest
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9.0];
julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];
julia> W = [2 1; 2 1; 2 1; 2 1; 2 1.0];
julia> P = [3 2; 3 2; 3 2; 3 2; 3 2.0];
julia> profit = deaprofit(X, Y, W, P, Gx = :Ones, Gy = :Ones)
julia> normfactor(profit)
5-element Vector{Float64}:
 8.0
 8.0
 8.0
 8.0
 8.0
```
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

# Examples
```jldoctest
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9.0];
julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];
julia> W = [2 1; 2 1; 2 1; 2 1; 2 1.0];
julia> P = [3 2; 3 2; 3 2; 3 2; 3 2.0];
julia> profit = deaprofit(X, Y, W, P, Gx = :Ones, Gy = :Ones, monetary = true);
julia> ismonetary(profit)
true
```
"""
function ismonetary(model::AbstractEconomicDEAModel)::Bool
    if isdefined(model, :monetary)
        return model.monetary
    else
        throw(ArgumentError("$(typeof(model)) has no monetary identifier"));
    end
end
