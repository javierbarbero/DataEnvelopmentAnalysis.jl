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
    dmunames::Union{Vector{AbstractString},Nothing}
    eff::Vector
    lambda::SparseMatrixCSC{Float64, Int64}
    crseff::Vector
    vrseff::Vector
    scaleff::Vector
    alloceff::Vector
    Xtarget::Matrix
    Ytarget::Matrix
end


"""
    deaprofitability(X, Y, W, P)
Compute profitability efficiency using data envelopment analysis for
inputs `X`, outputs `Y`, price of inputs `W`, and price of outputs `P`.

# Optional Arguments
- `alpha=0.5`: alpha to use for the generalized distance function.
- `names`: a vector of strings with the names of the decision making units.
"""
function deaprofitability(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector},
    W::Union{Matrix,Vector}, P::Union{Matrix,Vector};
    alpha::Float64 = 0.5,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::ProfitabilityDEAModel

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    nw, mw = size(W, 1), size(W, 2)
    np, sp = size(P, 1), size(P, 2)

    if nx != ny
        throw(DimensionMismatch("number of rows in X and Y ($nx, $ny) are not equal"));
    end
    if nw != nx
        throw(DimensionMismatch("number of rows in W and X ($nw, $nx) are not equal"));
    end
    if np != ny
        throw(DimensionMismatch("number of rows in P and Y ($np, $ny) are not equal"));
    end
    if mw != m
        throw(DimensionMismatch("number of columns in W and X ($mw, $m) are not equal"));
    end
    if sp != s
        throw(DimensionMismatch("number of columns in P and Y ($sp, $s) are not equal"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:NLP)
    end

    # Compute efficiency for each DMU
    n = nx

    Xtarget = zeros(n,m)
    Ytarget = zeros(n,s)
    pefficiency = zeros(n)
    plambdaeff = spzeros(n, n)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]
        w0 = W[i,:]
        p0 = P[i,:]

        # Create the optimization model
        deamodel = newdeamodel(optimizer)

        @variable(deamodel, eff, start = 1.0)
        @variable(deamodel, lambda[1:n] >= 0)

        @NLobjective(deamodel, Min, eff)

        @NLconstraint(deamodel, sum(sum(w0[mi] * X[t,mi] for mi in 1:m) / sum(p0[si] * Y[t,si] for si in 1:s) * lambda[t] for t in 1:n) == eff * sum(w0[j] * x0[j] for j in 1:m ) / sum(p0[j] * y0[j] for j in 1:s))

        @constraint(deamodel, sum(lambda) == 1)

        # Optimize and return results
        JuMP.optimize!(deamodel)

        pefficiency[i]  = JuMP.objective_value(deamodel)
        plambdaeff[i,:] = JuMP.value.(lambda)
        Xtarget[i,:] = X[i,:] .* pefficiency[i] ^(1-alpha)
        Ytarget[i,:] = Y[i,:] ./ ( pefficiency[i] ^alpha )

        # Check termination status
        if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    # Technical, scale and allocative efficiency
    crsefficiency = efficiency(deagdf(X, Y, alpha = alpha, rts = :CRS, slack = false, optimizer = optimizer))
    vrsefficiency = efficiency(deagdf(X, Y, alpha = alpha, rts = :VRS, slack = false, optimizer = optimizer))
    scalefficiency = crsefficiency ./ vrsefficiency
    allocefficiency = pefficiency ./ crsefficiency

    return ProfitabilityDEAModel(n, m, s, alpha, names, pefficiency, plambdaeff, crsefficiency, vrsefficiency, scalefficiency, allocefficiency, Xtarget, Ytarget)

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
        show(io, MIME"text/plain"(), CoefTable(hcat(eff, crseff, vrseff, scaleeff, alloceff), ["Profitability", "CRS", "VRS", "Scale", "Allocative"], dmunames))
    end

end
