```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Directional Distance Function Models

## Directional Distance Function Model

*Chambers, Chung and Fare (1996)* introduced a measure of efficiency that projects observation $\left( {{\mathbf{x}_o,\mathbf{y}_{o}}} \right)$
in a pre-assigned  direction  $\mathbf{g}= {\left({-{\mathbf{g_{x}^-},\mathbf{g^{+}_y}}} \right)\neq\mathbf{0}_{m+s}}$, $\mathbf{g^{-}_{x}}\mathbb{\in R}^m$ and  $\mathbf{g^{+}_{y}}\mathbb{\in R}^s$, in a proportion $\beta$. The associated linear program is:

```math
\begin{align}
\label{eq:ddf}
  & \underset{\beta ,\mathbf{\lambda }}{\mathop{\max }}\,\quad \quad \quad \quad \beta  \\
 & \text{subject}\ \text{to} \nonumber\\
 & \quad \quad \quad \quad \quad \ X\lambda\le {{\mathbf{x}}_{o}} -\beta{{\mathbf{g^-_x}}} \nonumber\\
 & \quad \quad \quad \quad \quad \  Y\mathbf{\lambda }\ge\ {{\mathbf{y}}_{o}}+\beta {{\mathbf{g^+_y}}}   \nonumber\\
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}.\nonumber\\  & \quad \nonumber
\end{align}
```

The measurement of technical efficiency assuming variable returns to scale, **VRS**, adds the following condition:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```

In this example we compute the directional distance function DEA model under constant returns to scale using ones as directions for both inputs and outputs:
```jldoctest 1
julia> using DataEnvelopmentAnalysis

julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaddf(X, Y, ones(size(X)), ones(size(Y)))
Directional DF DEA Model 
DMUs = 11; Inputs = 2; Outputs = 1
Returns to Scale = CRS
─────────────────────────────────────────────────────
      efficiency       slackX1       slackX2  slackY1
─────────────────────────────────────────────────────
1   -3.43053e-16   0.0           0.0              0.0
2    3.21996      -3.21359e-15   0.0              0.0
3    2.12169       0.0          -4.80367e-15      0.0
4    0.0          -8.03397e-16   0.0              0.0
5    6.73567      -2.41019e-15   0.0              0.0
6    1.94595      10.9189        0.0              0.0
7    0.0           0.0           0.0              0.0
8    3.63586       6.42718e-15   0.0              0.0
9    1.83784       4.75676       0.0              0.0
10  10.2311        6.12173e-15   0.0              0.0
11   0.0           0.0           4.0              0.0
─────────────────────────────────────────────────────
```

To compute the variable returns to scale model, we simply set the `rts` parameter to `:VRS`:
```jldoctest 1
julia> deaddf(X, Y, ones(size(X)), ones(size(Y)), rts = :VRS)
Directional DF DEA Model 
DMUs = 11; Inputs = 2; Outputs = 1
Returns to Scale = VRS
────────────────────────────────────────────────────
      efficiency       slackX1  slackX2      slackY1
────────────────────────────────────────────────────
1   -3.43053e-16   0.0              0.0  0.0        
2    1.41887       0.0              0.0  7.41268e-15
3    0.0           0.0              0.0  0.0        
4    0.0          -8.03397e-16      0.0  0.0        
5    4.06792       0.0              0.0  0.0        
6   -1.81673e-16   2.70127e-16      0.0  3.78178e-16
7    0.0           0.0              0.0  0.0        
8    0.0           0.0              0.0  0.0        
9    0.0           0.0              0.0  0.0        
10   5.0           0.0              6.0  0.0        
11   0.0           0.0              4.0  4.78849e-16
────────────────────────────────────────────────────
```

Estimated efficiency scores are returned with the `efficiency` function:
```jldoctest 1
julia> deaddfvrs = deaddf(X, Y, ones(size(X)), ones(size(Y)), rts = :VRS);

julia> efficiency(deaddfvrs)
11-element Array{Float64,1}:
 -3.4305304041327586e-16
  1.4188679245283022    
  0.0                   
  0.0                   
  4.067924528301886     
 -1.816728585750256e-16 
  0.0                   
  0.0                   
  0.0                   
  5.000000000000003     
  0.0     
```

The optimal peers, ``λ``, are returned with the `peers` function and are returned as a `SparseArrays.SparseMatrixCSC{Float64,Int64}` object:
```julia
peers(deaddfvrs)
```
### deaddf Function Documentation

```@docs
deaddf
```
