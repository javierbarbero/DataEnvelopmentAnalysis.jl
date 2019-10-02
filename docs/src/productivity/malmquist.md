```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# The Malmquist index

## The Malmquist Productivity Index

The Malmquist index introduced by *Caves, Christensen and Diewert(1982)* measures the change in  productivity of the observation under evaluation by comparing its relative performance with respect to  reference  technologies corresponding to two different time periods.

Following *Fare, Grosskopf, Norris and Zhang (1994)* productivity change can be decomposed into efficiency change and technical change under the assumption of a constant returns to scale techncology.

In this example we compute the Malmquist productivity index:
```jldoctest 1
julia> X = Array{Float64,3}(undef, 5, 1, 2);

julia> X[:, :, 1] = [2; 3; 5; 4; 4];

julia> X[:, :, 2] = [1; 2; 4; 3; 4];

julia> Y = Array{Float64,3}(undef, 5, 1, 2);

julia> Y[:, :, 1] = [1; 4; 6; 3; 5];

julia> Y[:, :, 2] = [1; 4; 6; 3; 3];

julia> malmquist(X, Y)
Mamlmquist DEA Model
DMUs = 5; Inputs = 1; Outputs = 1; Time periods = 2
Orientation = Input; Returns to Scale = CRS
Referene period = Geomean
─────────────────────────
         M        EC   TC
─────────────────────────
1  2.0      1.33333   1.5
2  1.5      1.0       1.5
3  1.25     0.833333  1.5
4  1.33333  0.888889  1.5
5  0.6      0.4       1.5
─────────────────────────
M  = Malmquist Productivity Index
EC = Efficiency Change
TC = Technological Change
```

### malmquist Function Documentation

```@docs
malmquist
```