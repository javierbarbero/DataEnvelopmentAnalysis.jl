# This file contains functions for the Cost Efficiency DEA model
"""
    CostDEAModel
An data structure representing a cost DEA model.
"""
struct CostDEAModel <: AbstractCostDEAModel
    n::Int64
    m::Int64
    s::Int64
    rts::Symbol
    disposY::Symbol
    dmunames::Union{Vector{String},Nothing}
    eff::Vector
    lambda::SparseMatrixCSC{Float64, Int64}
    techeff::Vector
    alloceff::Vector
    Xtarget::Matrix
    Ytarget::Matrix
end


"""
    deacost(X, Y, W)
Compute cost efficiency using data envelopment analysis for
inputs `X`, outputs `Y` and price of inputs `W`.

# Optional Arguments
- `rts=:VRS`: chooses variable returns to scale. For constant returns to scale choose `:CRS`.
- `dispos=:Strong`: chooses strong disposability of outputs. For weak disposability choose `:Weak`.
- `names`: a vector of strings with the names of the decision making units.

# Examples
```jldoctest
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9.0];

julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];

julia> W = [2 1; 2 1; 2 1; 2 1; 2 1.0];

julia> deacost(X, Y, W)
Cost DEA Model
DMUs = 5; Inputs = 2; Outputs = 2
Orientation = Input; Returns to Scale = VRS
──────────────────────────────────
       Cost  Technical  Allocative
──────────────────────────────────
1  0.615385      0.75     0.820513
2  1.0           1.0      1.0
3  1.0           1.0      1.0
4  0.5           0.5      1.0
5  0.347826      0.375    0.927536
──────────────────────────────────
```
"""
function deacost(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector},
    W::Union{Matrix,Vector}; rts::Symbol = :VRS, dispos::Symbol = :Strong,
    names::Union{Vector{String},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::CostDEAModel

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    nw, mw = size(W, 1), size(W, 2)

    if nx != ny
        throw(DimensionMismatch("number of rows in X and Y ($nx, $ny) are not equal"));
    end
    if nw != nx
        throw(DimensionMismatch("number of rows in W and X ($nw, $nx) are not equal"));
    end
    if mw != m
        throw(DimensionMismatch("number of columns in W and X ($mw, $m) are not equal"));
    end

    if dispos != :Strong && dispos != :Weak
        throw(ArgumentError("`disposX` must be :Strong or :Weak"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(GLPK.Optimizer)
    end

    # Get minimum cost targets and lambdas
    n = nx

    Xtarget, clambdaeff = deamincost(X, Y, W, rts = rts, dispos = dispos, optimizer = optimizer)
    Ytarget = Y[:,:]

    # Cost, technical and allocative efficiency
    cefficiency  = vec( sum(W .* Xtarget, dims = 2) ./ sum(W .* X, dims = 2) )
    techefficiency = efficiency(dea(X, Y, orient = :Input, rts = rts, slack = false, disposY = dispos, optimizer = optimizer))
    allocefficiency = cefficiency ./ techefficiency

    return CostDEAModel(n, m, s, rts, dispos, names, cefficiency, clambdaeff, techefficiency, allocefficiency, Xtarget, Ytarget)

end

ismonetary(model::CostDEAModel)::Bool = false;

function Base.show(io::IO, x::CostDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    disposY = x.disposY
    dmunames = names(x)

    eff = efficiency(x)
    techeff = efficiency(x, :Technical)
    alloceff = efficiency(x, :Allocative)

    if !compact
        print(io, "Cost DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = Input")
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        if disposY == :Weak print(io, "Weak disposability of outputs \n") end

        show(io, CoefTable(hcat(eff, techeff, alloceff), ["Cost", "Technical", "Allocative"], dmunames))
    end

end

