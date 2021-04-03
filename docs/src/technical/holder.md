```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Hölder Distance Function Models

*Briec (1998)* defined technical inefficiency using Hölder norms.

## Hölder 1

In this example we compute the Hölder l 1 DEA model under varible returns to scale:
```jldoctest 1
julia> using DataEnvelopmentAnalysis

julia> X = [2; 4; 8; 12; 6; 14; 14; 9.412];

julia> Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

julia> deaholder(X, Y, l = 1, rts = :VRS)
Hölder L1 DEA Model 
DMUs = 8; Inputs = 1; Outputs = 1
Orientation = Graph; Returns to Scale = VRS
────────────────────────────────────────────────
   efficiency  minimum      slackX1      slackY1
────────────────────────────────────────────────
1         0.0       X1  0.0          0.0
2         0.0       X1  0.0          0.0
3         0.0       X1  0.0          0.0
4         0.0       X1  0.0          0.0
5         3.0       X1  0.0          1.01506e-15
6         2.0       Y1  2.0          0.0
7         0.0       Y1  2.0          0.0
8         6.0       Y1  1.77636e-15  0.0
────────────────────────────────────────────────
```

Estimated efficiency scores are returned with the `efficiency` function:
```jldoctest 1
julia> holderl1 = deaholder(X, Y, l = 1, rts = :VRS);

julia> efficiency(holderl1)
8-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
 2.9999999999999996
 2.000000000000001
 0.0
 6.0
```

The input or output that determines the projection to the frontier is returned with:
```jldoctest 1
julia> efficiency(holderl1, :min)
8-element Vector{Int64}:
 1
 1
 1
 1
 1
 2
 2
 2
```
with inputs and outputs numbered sequentially.

## Hölder 2

!!! warning "Requieres a solver that supports SOS constraints"
    The Hölder l 2 model requieres a solver that supports SOS constraints, such as [Gurobi](https://github.com/jump-dev/Gurobi.jl). 

    Solving the model with Ipopt will return invalid results.

## Hölder Inf

In this example we compute the Hölder l Inf DEA model under varible returns to scale:
```jldoctest 1
julia> X = [2; 4; 8; 12; 6; 14; 14; 9.412];

julia> Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

julia> deaholder(X, Y, l = Inf, rts = :VRS)
Hölder LInf DEA Model 
DMUs = 8; Inputs = 1; Outputs = 1
Orientation = Graph; Returns to Scale = VRS
────────────────────────────────────────
   efficiency      slackX1       slackY1
────────────────────────────────────────
1       0.0    0.0           0.0
2       0.0    0.0           0.0
3       0.0    0.0           0.0
4       0.0    0.0           0.0
5       2.0    0.0          -8.08877e-16
6       2.0    0.0           0.0
7       0.0    2.0           0.0
8       3.832  2.96059e-16   0.0
────────────────────────────────────────
```

### deaholder Function Documentation

```@docs
deaholder
```
