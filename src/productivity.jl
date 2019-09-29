# This file contains types and structures for productivity DEA model
"""
    AbstractProductivityDEAlModel
An abstract type representing a productivity DEA model.
"""
abstract type AbstractProductivityDEAModel end

"""
    nbos(model::AbstractProductivityDEAModel)
Return number of observations of a productivity DEA model.
# Examples
```jldoctest

```
"""
nobs(model::AbstractProductivityDEAModel) = model.n

"""
    ninputs(model::AbstractProductivityDEAModel)
Return number of inputs of a productivity DEA model.
## Examples
```jldoctest

```
"""
ninputs(model::AbstractProductivityDEAModel) = model.m

"""
    noutputs(model::AbstractProductivityDEAModel)
Return number of outputs of a productivity DEA model.
# Examples
```jldoctest

```
"""
noutputs(model::AbstractProductivityDEAModel) = model.s

"""
    nperiods(model::AbstractProductivityDEAModel)
Return number of time periods of a productivity DEA model.
# Examples
```jldoctest

```
"""
nperiods(model::AbstractProductivityDEAModel) = model.periods

"""
    prodchange(model::AbstractProductivityDEAModel)
Return productivity change of a productivity change DEA model.

# Optional Arguments
- `type=Prod`: component of productivity change to return.

Type specification:
- `:Prod`: productivity change.
- `:EC`: efficiency change.
- `:TC`: technological change.

# Examples
```jldoctest
julia> X = Array{Float64,3}(undef, 5, 1, 2);

julia> X[:, :, 1] = [2; 3; 5; 4; 4];

julia> X[:, :, 2] = [1; 2; 4; 3; 4];

julia> Y = Array{Float64,3}(undef, 5, 1, 2);

julia> Y[:, :, 1] = [1; 4; 6; 3; 5];

julia> Y[:, :, 2] = [1; 4; 6; 3; 3];

julia> mprod = malmquist(X, Y);

julia> prodchange(mprod)
5×1 Array{Float64,2}:
 2.0
 1.5
 1.25
 1.3333333333333333
 0.6000000000000001

julia> prodchange(mprod, :EC)
5×1 Array{Float64,2}:
 1.3333333333333333
 1.0
 0.8333333333333334
 0.8888888888888888
 0.4

julia> prodchange(mprod, :TC)
5×1 Array{Float64,2}:
 1.5
 1.5
 1.5
 1.5
 1.5

```
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

    error(typeof(model), " has no productivity change component ", string(type))

end
