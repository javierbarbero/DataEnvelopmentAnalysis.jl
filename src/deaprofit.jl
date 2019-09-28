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
Compute profit efficinecy data envelopment analysis model for
inputs `X`, outputs `Y`, price of inputs `W`, and price of outptus `P`.

# Optional Arguments
- `Xref=X`: reference set of inputs to which evaluate the units.
- `Yref=Y`: reference set of outputs to which evaluate the units.

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
function deaprofit(X::Matrix, Y::Matrix, W::Matrix, P::Matrix, Gx::Matrix, Gy::Matrix; Xref::Matrix = X, Yref::Matrix = Y, Wref::Matrix = W, Pref::Matrix = P)::ProfitDEAModel
    # Check parameters
    nx, m = size(X)
    ny, s = size(Y)

    np, sp = size(P)

    nrefx, mref = size(Xref)
    nrefy, sref = size(Yref)

    if nx != ny
        error("number of observations is different in inputs and outputs")
    end
    if np != ny
        error("number of observations is different in output prices and outputs")
    end
    if sp != s
        error("number  of output prices and outputs is different")
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

    # Compute efficiency for each DMU
    n = nx
    nref = nrefx

    Xefficient = zeros(n,m)
    Yefficient = zeros(n,m)
    pefficiency = zeros(n)
    plambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        w0 = W[i,:]
        p0 = P[i,:]

        # Create the optimization model
        deamodel = Model(with_optimizer(GLPK.Optimizer))
        @variable(deamodel, Xeff[1:m])
        @variable(deamodel, Yeff[1:m])
        @variable(deamodel, lambda[1:nref] >= 0)

        @objective(deamodel, Max, (sum(p0[j] .* Yeff[j] for j in 1:s)) - (sum(w0[j] .* Xeff[j] for j in 1:m)))

        @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= Xeff[j])
        @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= Yeff[j])

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

    techefficiency = efficiency(deaddf(X, Y, Gx, Gy, rts = :VRS, Xref = Xref, Yref = Yref, slack = false))
    allocefficiency = pefficiency - techefficiency

    return ProfitDEAModel(n, m, s, pefficiency, plambdaeff, techefficiency, allocefficiency)

end

function deaprofit(X::Vector, Y::Matrix, W::Vector, P::Matrix, Gx::Vector, Gy::Matrix; Xref::Vector = X, Yref::Matrix = Y, Wref::Vector = W, Pref::Matrix = P)::ProfitDEAModel
    X = X[:,:]
    Xref = X[:,:]
    W = W[:,:]
    Wref = Wref[:,:]
    Gx = Gx[:,:]
    return deaprofit(X, Y, W, P, Gx, Gy, Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref)
end

function deaprofit(X::Matrix, Y::Vector, W::Matrix, P::Vector, Gx::Matrix, Gy::Vector; Xref::Matrix = X, Yref::Vector = Y, Wref::Matrix = W, Pref::Vector = P)::ProfitDEAModel
    Y = Y[:,:]
    Yref = Y[:,:]
    P = P[:,:]
    Pref = Pref[:,:]
    Gy = Gy[:,:]
    return deaprofit(X, Y, W, P, Gx, Gy, Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref)
end

function deaprofit(X::Vector, Y::Vector, W::Vector, P::Vector, Gx::Vector, Gy::Vector; Xref::Vector = X, Yref::Vector = Y, Wref::Vector = W, Pref::Vector = W)::ProfitDEAModel
    X = X[:,:]
    Xref = X[:,:]
    Y = Y[:,:]
    Yref = Y[:,:]
    W = W[:,:]
    Wref = Wref[:,:]
    P = P[:,:]
    Pref = Pref[:,:]
    Gx = Gx[:,:]
    Gy = Gy[:,:]
    return deaprofit(X, Y, W, P, Gx, Gy, Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref)
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
