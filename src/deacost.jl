# This file contains functions for the Cost Efficiency DEA model
"""
    CostDEAModel
An data structure representing a cost DEA model.
"""
struct CostDEAModel <: AbstractEconomicDEAModel
    n::Int64
    m::Int64
    s::Int64
    rts::Symbol
    eff::Vector
    lambda::SparseMatrixCSC{Float64, Int64}
    techeff::Vector
    alloceff::Vector
end


"""
    deacost(X, Y, W)
Compute cost efficiency using data envelopment analysis for
inputs `X`, outputs `Y` and price of inputs `W`.

# Optional Arguments
- `rts=:CRS`: chooses variable returns to scale. For constant returns to scale choose `:CRS`.

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
function deacost(X::Matrix, Y::Matrix, W::Matrix; rts::Symbol = :VRS, Xref::Matrix = X, Yref::Matrix = Y, Wref::Matrix = W)::CostDEAModel
    # Check parameters
    nx, m = size(X)
    ny, s = size(Y)

    nw, mw = size(W)

    nrefx, mref = size(Xref)
    nrefy, sref = size(Yref)

    if nx != ny
        error("number of observations is different in inputs and outputs")
    end
    if nw != nx
        error("number of observations is different in input prices and inputs")
    end
    if mw != m
        error("number  of input prices and intputs is different")
    end
    if nrefx != nrefy
        error("number of observations is different in inputs reference set and ouputs reference set")
    end
    if m != mref
        error("number of inputs in evaluation set and reference set is different")
    end
    if s != sref
        error("number of outputs in evaluation set and reference set is different")
    end
    if size(Wref) != size(Xref)
        error("size of reference prices for inputs should be equal to size of reference inputs")
    end

    # Compute efficiency for each DMU
    n = nx
    nref = nrefx

    Xefficient = zeros(n,m)
    cefficiency = zeros(n)
    clambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        y0 = Y[i,:]
        w0 = W[i,:]

        # Create the optimization model
        deamodel = Model(with_optimizer(GLPK.Optimizer))
        @variable(deamodel, Xeff[1:m])
        @variable(deamodel, lambda[1:nref] >= 0)

        @objective(deamodel, Min, sum(w0[j] .* Xeff[j] for j in 1:m))

        @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= Xeff[j])
        @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= y0[j])

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

        Xefficient[i,:]  = JuMP.value.(Xeff)
        clambdaeff[i,:] = JuMP.value.(lambda)

    end

    # Cost, technical and allocative efficiency
    cefficiency  = vec( sum(W .* Xefficient, dims = 2) ./ sum(W .* X, dims = 2) )
    techefficiency = efficiency(dea(X, Y, orient = :Input, rts = rts, Xref = Xref, Yref = Yref, slack = false))
    allocefficiency = cefficiency ./ techefficiency
    return CostDEAModel(n, m, s, rts, cefficiency, clambdaeff, techefficiency, allocefficiency)

end

function deacost(X::Vector, Y::Matrix, W::Vector, rts::Symbol = :VRS, Xref::Vector = X, Yref::Matrix = Y, Wref::Vector = W)::CostDEAModel
    X = X[:,:]
    Xref = Xref[:,:]
    W = W[:,:]
    Wref = Wref[:,:]
    return deacost(X, Y, W, rts = rts, Xref = Xref, Yref = Yref, Wref = Wref)
end

function deacost(X::Matrix, Y::Vector, W::Matrix; rts::Symbol = :VRS, Xref::Matrix = X, Yref::Vector = Y, Wref::Matrix = W)::CostDEAModel
    Y = Y[:,:]
    Yref = Yref[:,:]
    return deacost(X, Y, W, rts = rts, Xref = Xref, Yref = Yref, Wref = Wref)
end

function deacost(X::Vector, Y::Vector, W::Vector; rts::Symbol = :VRS, Xref::Vector = X, Yref::Vector = Y, Wref::Vector = W)::CostDEAModel
    X = X[:,:]
    Xref = Xref[:,:]
    W = W[:,:]
    Wref = Wref[:,:]
    Y = Y[:,:]
    Yref = Yref[:,:]
    return deacost(X, Y, W, rts = rts, Xref = Xref, Yref = Yref, Wref = Wref)
end

function Base.show(io::IO, x::CostDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
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
        show(io, CoefTable(hcat(eff, techeff, alloceff), ["Cost", "Technical", "Allocative"], ["$i" for i in 1:n]))

    else

    end
end
