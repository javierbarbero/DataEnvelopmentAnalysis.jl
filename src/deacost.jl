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
    names::Union{Vector{String},Nothing} = nothing)::CostDEAModel

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    nw, mw = size(W, 1), size(W, 2)

    if nx != ny
        error("number of observations is different in inputs and outputs")
    end
    if nw != nx
        error("number of observations is different in input prices and inputs")
    end
    if mw != m
        error("number of input prices and intputs is different")
    end

    if dispos != :Strong && dispos != :Weak
        error("Invalued disposability $dispos. Disposability should be :Strong or :Weak")
    end

    # Compute efficiency for each DMU
    n = nx

    Xtarget = zeros(n,m)
    Ytarget = Y[:,:]
    cefficiency = zeros(n)
    clambdaeff = spzeros(n, n)

    for i=1:n
        # Value of inputs and outputs to evaluate
        y0 = Y[i,:]
        w0 = W[i,:]

        # Create the optimization model
        deamodel = Model(GLPK.Optimizer)

        @variable(deamodel, Xeff[1:m])
        @variable(deamodel, lambda[1:n] >= 0)

        @objective(deamodel, Min, sum(w0[j] .* Xeff[j] for j in 1:m))

        @constraint(deamodel, [j in 1:m], sum(X[t,j] * lambda[t] for t in 1:n) <= Xeff[j])

        if dispos == :Strong
            @constraint(deamodel, [j in 1:s], sum(Y[t,j] * lambda[t] for t in 1:n) >= y0[j])
        elseif dispos == :Weak
            @constraint(deamodel, [j in 1:s], sum(Y[t,j] * lambda[t] for t in 1:n) == y0[j])
        end

        # Add return to scale constraints
        if rts == :CRS
            # No contraint to add for constant returns to scale
        elseif rts == :VRS
            @constraint(deamodel, sum(lambda) == 1)
        else
            error("Invalid returns to scale $rts. Returns to scale should be :CRS or :VRS")
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        Xtarget[i,:]  = JuMP.value.(Xeff)
        clambdaeff[i,:] = JuMP.value.(lambda)

        # Check termination status
        if termination_status(deamodel) != MOI.OPTIMAL
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    # Cost, technical and allocative efficiency
    cefficiency  = vec( sum(W .* Xtarget, dims = 2) ./ sum(W .* X, dims = 2) )
    techefficiency = efficiency(dea(X, Y, orient = :Input, rts = rts, slack = false, disposY = dispos))
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

