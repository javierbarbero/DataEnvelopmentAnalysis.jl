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

