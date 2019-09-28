# This file contains types and structures for technical DEA model
"""
    AbstractTechnicaDEAlModel
An abstract type representing a technical DEA model.
"""
abstract type AbstractTechnicalDEAModel end

"""
    efficiency(model::AbstractTechnicalDEAModel)
Return efficiency scores of a technical DEA model.
# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> efficiency(deaio)
11-element Array{Float64,1}:
 1.0
 0.6222896790980051
 0.8198562443845464
 1.0
 0.3103709311127934
 0.5555555555555556
 1.0
 0.7576690895651103
 0.8201058201058201
 0.49056603773584917
 1.0
```
"""
efficiency(model::AbstractTechnicalDEAModel) = model.eff

"""
    nbos(model::AbstractTechnicalDEAModel)
Return number of observations of a technical DEA model.
# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> nobs(deaio)
11
```
"""
nobs(model::AbstractTechnicalDEAModel) = model.n

"""
    ninputs(model::AbstractTechnicalDEAModel)
Return number of inputs of a technical DEA model.
## Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> ninputs(deaio)
2
```
"""
ninputs(model::AbstractTechnicalDEAModel) = model.m

"""
    noutputs(model::AbstractTechnicalDEAModel)
Return number of outputs of a technical DEA model.
# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> noutputs(deaio)
1
```
"""
noutputs(model::AbstractTechnicalDEAModel) = model.s

"""
    peers(model::AbstractTechnicalDEAModel)
Return peers of a technical DEA model.
# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> peers(deaio)
11×11 SparseArrays.SparseMatrixCSC{Float64,Int64} with 17 stored entries:
  [1 ,  1]  =  1.0
  [3 ,  1]  =  1.13432
  [2 ,  4]  =  0.424978
  [3 ,  4]  =  0.438005
  [4 ,  4]  =  1.0
  ⋮
  [7 ,  7]  =  1.0
  [8 ,  7]  =  0.114574
  [9 ,  7]  =  1.14815
  [10,  7]  =  0.490566
  [11, 11]  =  1.0
```
"""
peers(model::AbstractTechnicalDEAModel) = model.lambda

"""
    slacks(model::AbstractTechnicalDEAModel, slack::Symbol)
Return slacks a technical DEA model.
# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> slacks(deaio, :X)
11×2 Array{Float64,2}:
  0.0          0.0        
 -4.41868e-15  0.0        
  0.0          8.17926e-15
 -8.03397e-16  0.0        
  1.80764e-15  0.0        
  4.44444      0.0        
  0.0          0.0        
  1.60679e-15  0.0        
  1.64021      0.0        
  9.68683e-15  0.0        
  0.0          4.0        

julia> slacks(deaio, :Y)
11×1 Array{Float64,2}:
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
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

end