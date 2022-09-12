```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Configuring the optimizer

DataEnvelopmentAnalysis.jl will use a default optimizer/solver for each DEA model, as shown in the next table.

| Function       | Specific Options | Problem type | Default Optimizer |
| ---------------|--------:|------------:|------------------:| 
| `dea`          |         | LP           | GLPK              |
| `deam`         |         | LP           | GLPK              |
| `deaboot`      |         | LP           | GLPK              |
| `deabigdata`   |         | LP           | GLPK              |
| `deaddf`       |         | LP           | GLPK              |
| `deaadd`       |         | LP           | GLPK              |
| `deagdf`       |         | NLP          | Ipopt             |
| `dearussell`   | `:Input` or `:Output`        | LP           | GLPK              |
| `dearussell`   | `:Graph`        | NLP     | Ipopt      |
| `deaerg`       |         | LP           | GLPK              |
| `deamddf`      |         | LP           | GLPK              |
| `deaholder`    | `l = 1` | LP           | GLPK              |
| `deaholder`    | `l = 2` | QP           |                   |
| `deaholder`    |`l = Inf`| LP           | GLPK              |
| `dearddf`      | `:ERG`  | LP           | GLPK              |
| `dearddf`      | `:MDDF` | LP           | GLPK              |
| `deacost`      |         | LP           | GLPK              |
| `dearevenue`   |         | LP           | GLPK              |
| `deaprofit`    |         | LP           | GLPK              |
| `deaprofitability` |         | NLP          | Ipopt         |
| `malmquist`    |          |LP           | GLPK              |


Where:
- LP = Linear programming.
- NLP = Nonlinear programming.
- QP = Quadratic programming.

Models can be solved using a different optimizer by passing a `DEAOptimizer` object to the `optimizer` optional argument. See [JuMP documentation](https://jump.dev/JuMP.jl/v0.21.6/installation/#Installing-a-solver) for a list of all available solvers.

!!! warning "Choose a valid optimizer"
    The optimizer must support the problem type of the DEA model.

    For example, you cannot solve a Generalized Distance Function DEA model using the GLPK solver because it is a linear programming solver and `deagdf` requires a nonlinear programming solver.

The following is an example of solving the radial DEA model using the `Ipopt` sovler:
```jldoctest
julia> using Ipopt

julia> using DataEnvelopmentAnalysis

julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> myoptimizer = DEAOptimizer(Ipopt.Optimizer, time_limit = 10, silent = true);

julia> dea(X, Y, slack = false, optimizer = myoptimizer)
Radial DEA Model 
DMUs = 11; Inputs = 2; Outputs = 1
Orientation = Input; Returns to Scale = CRS
──────────────
    efficiency
──────────────
1     1.0
2     0.62229
3     0.819856
4     1.0
5     0.310371
6     0.555555
7     1.0
8     0.757669
9     0.820106
10    0.490566
11    1.0
──────────────
```

### Optimizer API

```@docs
DEAOptimizer
newdeamodel
```
