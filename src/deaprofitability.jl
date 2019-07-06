# This file contains functions for the Profitability Efficiency DEA model
"""
    ProfitabilityDEAModel
An data structure representing a profitability DEA model.
"""
struct ProfitabilityDEAModel <: AbstractEconomicDEAModel
    n::Int64
    m::Int64
    s::Int64
    alpha::Float64
    eff::Vector
    lambda::SparseMatrixCSC{Float64, Int64}
    crseff::Vector
    vrseff::Vector
    scaleff::Vector
    alloceff::Vector
end


"""
    deaprofitability(X, Y, W, P)
Computes profitability efficinecy data envelopment analysis model for
inputs `X`, outputs `Y`, price of inputs `W`, and price of outptus `P`.

### Optional Arguments
- `alpha=0.5`: alpha to use for the generalized distance function.
- `Xref=X`: reference set of inputs to which evaluate the units.
- `Yref=Y`: reference set of outputs to which evaluate the units.

# Examples
```jldoctest
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9.0];
julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];
julia> W = [2 1; 2 1; 2 1; 2 1; 2 1.0];
julia> P = [3 2; 3 2; 3 2; 3 2; 3 2.0];
julia> deaprofitability(X, Y, W, P)
Profitability DEA Model
DMUs = 5; Inputs = 2; Outputs = 2
alpha = 0.5; Returns to Scale = VRS
─────────────────────────────────────────────────────────
   Profitability       CRS      VRS     Scale  Allocative
─────────────────────────────────────────────────────────
1       0.38796   0.636364  0.68185  0.93329     0.609651
2       1.0       1.0       1.0      1.0         1.0
3       0.765217  1.0       1.0      1.0         0.765217
4       0.25      0.25      0.25     1.0         1.0
5       0.15879   0.26087   0.36     0.724638    0.608696
─────────────────────────────────────────────────────────
```
"""
function deaprofitability(X::Matrix, Y::Matrix, W::Matrix, P::Matrix; alpha::Float64 = 0.5, Xref::Matrix = X, Yref::Matrix = Y, Wref::Matrix = W, Pref::Matrix = P)::ProfitabilityDEAModel
    # Check parameters
    nx, m = size(X)
    ny, s = size(Y)

    nw, mw = size(W)
    np, mp = size(P)

    nrefx, mref = size(Xref)
    nrefy, sref = size(Yref)

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
        error("number  of input prices and intputs is different")
    end
    if mp != m
        error("number of output prices and otuputs is different")
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

    pefficiency = zeros(n)
    plambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]
        w0 = W[i,:]
        p0 = P[i,:]

        # Create the optimization model
        deamodel = Model(with_optimizer(Ipopt.Optimizer, print_level= 0 ))
        @variable(deamodel, eff, start = 1.0)
        @variable(deamodel, lambda[1:nref] >= 0)

        @NLobjective(deamodel, Min, eff)

        @NLconstraint(deamodel, sum(sum(Wref[t,mi] * Xref[t,mi] for mi in 1:m) / sum(Pref[t,si] * Yref[t,si] for si in 1:s) * lambda[t] for t in 1:nref) == eff * sum(w0[j] * x0[j] for j in 1:m ) / sum(p0[j] * y0[j] for j in 1:s))

        @constraint(deamodel, sum(lambda) == 1)

        # Optimize and return results
        JuMP.optimize!(deamodel)

        pefficiency[i]  = JuMP.objective_value(deamodel)
        plambdaeff[i,:] = JuMP.value.(lambda)

    end

    # Technical efficiency
    crsefficiency = efficiency(deagdf(X, Y, alpha, rts = :CRS, Xref = Xref, Yref = Yref))
    vrsefficiency = efficiency(deagdf(X, Y, alpha, rts = :VRS, Xref = Xref, Yref = Yref))
    scalefficiency = crsefficiency ./ vrsefficiency
    allocefficiency = pefficiency ./ crsefficiency

    return ProfitabilityDEAModel(n, m, s, alpha, pefficiency, plambdaeff, crsefficiency, vrsefficiency, scalefficiency, allocefficiency)

end

function deaprofitability(X::Vector, Y::Matrix, W::Vector, P::Matrix; alpha::Float64 = 0.5, Xref::Vector = X, Yref::Matrix = Y, Wref::Vector = W, Pref::Matrix = P)::ProfitabilityDEAModel
    X = X[:,:]
    Xref = X[:,:]
    W = W[:,:]
    Wref = Wref[:,:]
    return deaprofitability(X, Y, W, P, alpha = alpha, Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref)
end

function deaprofitability(X::Matrix, Y::Vector, W::Matrix, P::Vector; alpha::Float64 = 0.5, Xref::Matrix = X, Yref::Vector = Y, Wref::Matrix = W, Pref::Vector = P)::ProfitabilityDEAModel
    Y = Y[:,:]
    Yref = Y[:,:]
    P = P[:,:]
    Pref = Pref[:,:]
    return deaprofitability(X, Y, W, P, alpha = alpha, Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref)
end

function deaprofitability(X::Vector, Y::Vector, W::Vector, P::Vector; alpha::Float64 = 0.5, Xref::Vector = X, Yref::Vector = Y, Wref::Vector = W, Pref::Vector = P)::ProfitabilityDEAModel
    X = X[:,:]
    Xref = X[:,:]
    W = W[:,:]
    Wref = Wref[:,:]
    Y = Y[:,:]
    Yref = Y[:,:]
    P = P[:,:]
    Pref = Pref[:,:]
    return deaprofitability(X, Y, W, P, alpha = alpha, Xref = Xref, Yref = Yref, Wref = Wref, Pref = Pref)
end

function Base.show(io::IO, x::ProfitabilityDEAModel)
    compact = get(io, :compact, false)

    eff = efficiency(x)
    crseff = efficiency(x, :CRS)
    vrseff = efficiency(x, :VRS)
    scaleeff = efficiency(x, :Scale)
    alloceff = efficiency(x, :Allocative)
    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)

    if !compact
        print(io, "Profitability DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "alpha = ", x.alpha)
        print(io, "; Returns to Scale = VRS")
        print(io, "\n")
        show(io, CoefTable(hcat(eff, crseff, vrseff, scaleeff, alloceff), ["Profitability", "CRS", "VRS", "Scale", "Allocative"], ["$i" for i in 1:n]))
    else

    end
end
