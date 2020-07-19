# This file contains functions for the Profitability Efficiency DEA model
"""
    ProfitabilityDEAModel
An data structure representing a profitability DEA model.
"""
struct ProfitabilityDEAModel <: AbstractProfitabilityDEAModel
    n::Int64
    m::Int64
    s::Int64
    alpha::Float64
    dmunames::Union{Vector{String},Nothing}
    eff::Vector
    lambda::SparseMatrixCSC{Float64, Int64}
    crseff::Vector
    vrseff::Vector
    scaleff::Vector
    alloceff::Vector
end


"""
    deaprofitability(X, Y, W, P)
Compute profitability efficiency using data envelopment analysis for
inputs `X`, outputs `Y`, price of inputs `W`, and price of outputs `P`.

# Optional Arguments
- `alpha=0.5`: alpha to use for the generalized distance function.
- `names`: a vector of strings with the names of the decision making units.

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
function deaprofitability(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector},
    W::Union{Matrix,Vector}, P::Union{Matrix,Vector};
    alpha::Float64 = 0.5,
    names::Union{Vector{String},Nothing} = nothing)::ProfitabilityDEAModel

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    nw, mw = size(W, 1), size(W, 2)
    np, sp = size(P, 1), size(P, 2)

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
        error("number of output prices and otuputs is different")
    end

    # Compute efficiency for each DMU
    n = nx

    pefficiency = zeros(n)
    plambdaeff = spzeros(n, n)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]
        w0 = W[i,:]
        p0 = P[i,:]

        # Create the optimization model
        deamodel = Model(Ipopt.Optimizer)
        set_silent(deamodel)
        set_optimizer_attribute(deamodel, "print_level", 0)

        @variable(deamodel, eff, start = 1.0)
        @variable(deamodel, lambda[1:n] >= 0)

        @NLobjective(deamodel, Min, eff)

        @NLconstraint(deamodel, sum(sum(w0[mi] * X[t,mi] for mi in 1:m) / sum(p0[si] * Y[t,si] for si in 1:s) * lambda[t] for t in 1:n) == eff * sum(w0[j] * x0[j] for j in 1:m ) / sum(p0[j] * y0[j] for j in 1:s))

        @constraint(deamodel, sum(lambda) == 1)

        # Optimize and return results
        JuMP.optimize!(deamodel)

        pefficiency[i]  = JuMP.objective_value(deamodel)
        plambdaeff[i,:] = JuMP.value.(lambda)

        # Check termination status
        if termination_status(deamodel) != MOI.LOCALLY_SOLVED
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    # Technical, scale and allocative efficiency
    crsefficiency = efficiency(deagdf(X, Y, alpha = alpha, rts = :CRS, slack = false))
    vrsefficiency = efficiency(deagdf(X, Y, alpha = alpha, rts = :VRS, slack = false))
    scalefficiency = crsefficiency ./ vrsefficiency
    allocefficiency = pefficiency ./ crsefficiency

    return ProfitabilityDEAModel(n, m, s, alpha, names, pefficiency, plambdaeff, crsefficiency, vrsefficiency, scalefficiency, allocefficiency)

end

function Base.show(io::IO, x::ProfitabilityDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    dmunames = names(x)

    eff = efficiency(x)
    crseff = efficiency(x, :CRS)
    vrseff = efficiency(x, :VRS)
    scaleeff = efficiency(x, :Scale)
    alloceff = efficiency(x, :Allocative)

    if !compact
        print(io, "Profitability DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "alpha = ", x.alpha)
        print(io, "; Returns to Scale = VRS")
        print(io, "\n")
        show(io, CoefTable(hcat(eff, crseff, vrseff, scaleeff, alloceff), ["Profitability", "CRS", "VRS", "Scale", "Allocative"], dmunames))
    end

end
