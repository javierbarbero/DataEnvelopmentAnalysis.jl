# This file contains functions for the Revenue Efficiency DEA model
"""
    RevenueDEAModel
An data structure representing a revenue DEA model.
"""
struct RevenueDEAModel <: AbstractRevenueDEAModel
    n::Int64
    m::Int64
    s::Int64
    rts::Symbol
    disposX::Symbol
    dmunames::Union{Vector{String},Nothing}
    eff::Vector
    lambda::SparseMatrixCSC{Float64, Int64}
    techeff::Vector
    alloceff::Vector
    Xtarget::Matrix
    Ytarget::Matrix
end


"""
    dearevenue(X, Y, P)
Compute revenue efficiency using data envelopment analysis for
inputs `X`, outputs `Y` and price of outputs `P`.

# Optional Arguments
- `rts=:VRS`: chooses variable returns to scale. For constant returns to scale choose `:CRS`.
- `dispos=:Strong`: chooses strong disposability of inputs. For weak disposability choose `:Weak`.
- `names`: a vector of strings with the names of the decision making units.

# Examples
```jldoctest
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9.0];

julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];

julia> P = [3 2; 3 2; 3 2; 3 2; 3 2.0];

julia> dearevenue(X, Y, P)
Revenue DEA Model
DMUs = 5; Inputs = 2; Outputs = 2
Orientation = Output; Returns to Scale = VRS
──────────────────────────────────
    Revenue  Technical  Allocative
──────────────────────────────────
1  0.644444   0.777778    0.828571
2  1.0        1.0         1.0
3  1.0        1.0         1.0
4  0.5        0.5         1.0
5  0.456522   0.6         0.76087
──────────────────────────────────
```
"""
function dearevenue(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector},
    P::Union{Matrix,Vector}; rts::Symbol = :VRS, dispos::Symbol = :Strong,
    names::Union{Vector{String},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::RevenueDEAModel

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    np, sp = size(P, 1), size(P, 2)

    if nx != ny
        error("number of observations is different in inputs and outputs")
    end
    if np != ny
        error("number of observations is different in output prices and outputs")
    end
    if sp != s
        error("number of output prices and outputs is different")
    end

    if dispos != :Strong && dispos != :Weak
        error("Invalued disposability $dispos. Disposability should be :Strong or :Weak")
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(GLPK.Optimizer)
    end    

    # Get maximum revenue targets and lambdas
    n = nx

    Xtarget = X[:,:]
    Ytarget, rlambdaeff = deamaxrevenue(X, Y, P, rts = rts, dispos = dispos, optimizer = optimizer)

    # Revenue, technical and allocative efficiency
    refficiency  = vec( sum(P .* Y, dims = 2) ./ sum(P .* Ytarget, dims = 2) )
    techefficiency = 1 ./ efficiency(dea(X, Y, orient = :Output, rts = rts, slack = false, disposX = dispos, optimizer = optimizer))
    allocefficiency = refficiency ./ techefficiency
    return RevenueDEAModel(n, m, s, rts, dispos, names, refficiency, rlambdaeff, techefficiency, allocefficiency, Xtarget, Ytarget)

end

ismonetary(model::RevenueDEAModel)::Bool = false;

function Base.show(io::IO, x::RevenueDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    disposX = x.disposX
    dmunames = names(x)

    eff = efficiency(x)
    techeff = efficiency(x, :Technical)
    alloceff = efficiency(x, :Allocative)

    if !compact
        print(io, "Revenue DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = Output")
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        if disposX == :Weak print(io, "Weak disposability of inputs \n") end

        show(io, CoefTable(hcat(eff, techeff, alloceff), ["Revenue", "Technical", "Allocative"], dmunames))
    end

end
