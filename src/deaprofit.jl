# This file contains functions for the Profit Efficiency DEA model
"""
    ProfitDEAModel
An data structure representing a profit DEA model.
"""
struct ProfitDEAModel <: AbstractEconomicDEAModel
    n::Int64
    m::Int64
    s::Int64
    eff::Vector
    lambda::SparseMatrixCSC{Float64, Int64}
    techeff::Vector
    alloceff::Vector
end


"""
    deaprofit(X, Y, P)
Compute profit efficiency using data envelopment analysis model for
inputs `X`, outputs `Y`, price of inputs `W`, and price of outputs `P`.

# Examples
```jldoctest
julia> X = [1 1; 1 1; 0.75 1.5; 0.5 2; 0.5 2; 2 2; 2.75 3.5; 1.375 1.75];

julia> Y = [1 11; 5 3; 5 5; 2 9; 4 5; 4 2; 3 3; 4.5 3.5];

julia> P = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1];

julia> W = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1];

julia> GxGydollar = 1 ./ (sum(P, dims = 2) + sum(W, dims = 2));

julia> Gx = repeat(GxGydollar, 1, 2);

julia> Gy = repeat(GxGydollar, 1, 2);

julia> deaprofit(X, Y, W, P, Gx, Gy)
Profit DEA Model 
DMUs = 8; Inputs = 2; Outputs = 2
Returns to Scale = VRS
─────────────────────────────────────
   Profit     Technical    Allocative
─────────────────────────────────────
1     2.0   0.0           2.0        
2     2.0  -5.41234e-16   2.0        
3     0.0   0.0           0.0        
4     2.0   0.0           2.0        
5     2.0   0.0           2.0        
6     8.0   6.0           2.0        
7    12.0  12.0          -1.77636e-15
8     4.0   3.0           1.0        
─────────────────────────────────────
```
"""
function deaprofit(X::Matrix, Y::Matrix, W::Matrix, P::Matrix, Gx::Matrix, Gy::Matrix)::ProfitDEAModel
    # Check parameters
    nx, m = size(X)
    ny, s = size(Y)

    nw, mw = size(W)
    np, sp = size(P)

    if nx != ny
        error("number of observations is different in inputs and outputs")
    end
    if nw != nx
        error("number of observations is different in input prices and inputs")
    end
    if np != ny
        error("number of observations is different in output prices and outputs")
    end
    if mw != m
        error("number of input prices and intputs is different")
    end
    if sp != s
        error("number of output prices and outputs is different")
    end
    if size(Gx) != size(X)
        error("size of inputs should be equal to size of inputs direction")
    end
    if size(Gy) != size(Y)
        error("size of outputs should be equal to size of outputs direction")
    end

    # Compute efficiency for each DMU
    n = nx

    Xefficient = zeros(n,m)
    Yefficient = zeros(n,m)
    pefficiency = zeros(n)
    plambdaeff = spzeros(n, n)

    for i=1:n
        # Value of inputs and outputs to evaluate
        w0 = W[i,:]
        p0 = P[i,:]

        # Create the optimization model
        deamodel = Model(with_optimizer(GLPK.Optimizer))
        @variable(deamodel, Xeff[1:m])
        @variable(deamodel, Yeff[1:m])
        @variable(deamodel, lambda[1:n] >= 0)

        @objective(deamodel, Max, (sum(p0[j] .* Yeff[j] for j in 1:s)) - (sum(w0[j] .* Xeff[j] for j in 1:m)))

        @constraint(deamodel, [j in 1:m], sum(X[t,j] * lambda[t] for t in 1:n) <= Xeff[j])
        @constraint(deamodel, [j in 1:s], sum(Y[t,j] * lambda[t] for t in 1:n) >= Yeff[j])

        @constraint(deamodel, sum(lambda) == 1)

        # Optimize and return results
        JuMP.optimize!(deamodel)

        Xefficient[i,:]  = JuMP.value.(Xeff)
        Yefficient[i,:]  = JuMP.value.(Yeff)
        plambdaeff[i,:] = JuMP.value.(lambda)

    end

    # Profit, technical and allocative efficiency
    maxprofit = sum(P .* Yefficient, dims = 2) .- sum(W .* Xefficient, dims = 2)

    pefficiency_num  = maxprofit .- ( sum(P .* Y, dims = 2) .- sum(W .* X, dims = 2))
    pefficiency_den = sum(P .* Gy, dims = 2) .+ sum(W .* Gx, dims = 2)
    pefficiency = vec( pefficiency_num ./ pefficiency_den )

    techefficiency = efficiency(deaddf(X, Y, Gx, Gy, rts = :VRS, slack = false))
    allocefficiency = pefficiency - techefficiency

    return ProfitDEAModel(n, m, s, pefficiency, plambdaeff, techefficiency, allocefficiency)

end

function deaprofit(X::Vector, Y::Matrix, W::Vector, P::Matrix, Gx::Vector, Gy::Matrix)::ProfitDEAModel
    X = X[:,:]
    W = W[:,:]
    Gx = Gx[:,:]
    return deaprofit(X, Y, W, P, Gx, Gy)
end

function deaprofit(X::Matrix, Y::Vector, W::Matrix, P::Vector, Gx::Matrix, Gy::Vector)::ProfitDEAModel
    Y = Y[:,:]
    P = P[:,:]
    Gy = Gy[:,:]
    return deaprofit(X, Y, W, P, Gx, Gy)
end

function deaprofit(X::Vector, Y::Vector, W::Vector, P::Vector, Gx::Vector, Gy::Vector)::ProfitDEAModel
    X = X[:,:]
    Y = Y[:,:]
    W = W[:,:]
    P = P[:,:]
    Gx = Gx[:,:]
    Gy = Gy[:,:]
    return deaprofit(X, Y, W, P, Gx, Gy)
end

function Base.show(io::IO, x::ProfitDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    eff = efficiency(x)
    techeff = efficiency(x, :Technical)
    alloceff = efficiency(x, :Allocative)

    if !compact
        print(io, "Profit DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Returns to Scale = VRS")
        print(io, "\n")
        show(io, CoefTable(hcat(eff, techeff, alloceff), ["Profit", "Technical", "Allocative"], ["$i" for i in 1:n]))

    else

    end
end
