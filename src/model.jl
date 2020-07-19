# This file contains types and structures for DEA model
"""
    AbstractDEAModel
An abstract type representing a DEA model.
"""
abstract type AbstractDEAModel end

"""
    nbos(model::Abstract DEAModel)
Return number of observations of a DEA model.
# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> nobs(deaio)
11
```
"""
nobs(model::AbstractDEAModel) = model.n

"""
    ninputs(model::AbstractDEAModel)
Return number of inputs of a DEA model.
## Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> ninputs(deaio)
2
```
"""
ninputs(model::AbstractDEAModel) = model.m

"""
    noutputs(model::AbstractDEAModel)
Return number of outputs of a DEA model.
# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> noutputs(deaio)
1
```
"""
noutputs(model::AbstractDEAModel) = model.s

"""
    names(model::AbstractDEAModel)
Return the names of the decision making units (DMUs)
# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> dmunames = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]

julia> deaio = dea(X, Y, names = dmunames);

julia> names(deaio)
11-element Array{String,1}:
 "A"
 "B"
 "C"
 "D"
 "E"
 "F"
 "G"
 "H"
 "I"
 "J"
 "K"
```
"""
function Base.names(model::AbstractDEAModel)

    xnobs = nobs(model)

    if isdefined(model, :dmunames)

        if model.dmunames === nothing
            # If model have no names, return numeric sequence
            retnames = ["$i" for i in 1:xnobs]
        else
            xnameslength = length(model.dmunames)

            if xnameslength == xnobs
                retnames = model.dmunames
            elseif xnameslength < xnobs
                # If length of names is lower than number of observations, append numbers to match
                @warn("Length of names lower than number of observations")
                retnames = [model.dmunames; ["$i" for i in (xnameslength + 1):xnobs]]
            elseif xnameslength > xnobs
                # If length of names is greater than number of observations, split
                @warn("Length of names greater than number of observations")
                retnames = model.dmunames[1:xnobs]
            end
        end
    else
        # Return numeric sequence
        retnames = ["$i" for i in 1:xnobs]
    end

    return retnames
end
