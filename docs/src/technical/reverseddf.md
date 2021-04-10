```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Reverse Directional Distance Function

In this example, we compute the Reverse Directional Distance Function (Pastor et al., 2016) DEA model for the Enhanced Russell Graph associated efficiency measure under variable returns to scale:
```jldoctest 1
julia> X = [2; 4; 8; 12; 6; 14; 14; 9.412];

julia> Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

julia> dearddf(X, Y, :ERG, rts = :VRS)
Reverse DDF DEA Model 
DMUs = 8; Inputs = 1; Outputs = 1
Orientation = Graph; Returns to Scale = VRS
Associated efficiency measure = ERG
─────────────
   efficiency
─────────────
1    0.0
2    0.0
3    0.0
4    0.0
5    0.6
6    0.52381
7    0.142857
8    0.8
─────────────
```

Estimated efficiency scores are returned with the `efficiency` function:
```jldoctest 1
julia> rddferg = dearddf(X, Y, :ERG, rts = :VRS);

julia> efficiency(rddferg)
8-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
 0.6000000000000002
 0.523809523809524
 0.14285714285714296
 0.8
```

### dearddf Function Documentation

```@docs
dearddf
```
