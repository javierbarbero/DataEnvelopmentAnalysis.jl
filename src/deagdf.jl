# This file contains functions for the Generalized Distance Function DEA model
"""
    GeneralizedDFDEAModel
An data structure representing a generalized distance function DEA model.
"""
struct GeneralizedDFDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64
    alpha::Float64
    rts::Symbol
    dmunames::Union{Vector{String},Nothing}
    eff::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    Xtarget::Matrix
    Ytarget::Matrix
end

"""
    deagdf(X, Y, alpha)
Compute generalized distance function data envelopment analysis model for
inputs `X`, outputs `Y`, and `alpha`.

# Optional Arguments
- `alpha=0.5`: alpha value.
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `slack=true`: compute input and output slacks.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `names`: a vector of strings with the names of the decision making units.

# Examples
```jldoctest
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9];

julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6];

julia> deagdf(X, Y, alpha = 0.5, rts = :VRS)
Generalized DF DEA Model 
DMUs = 5; Inputs = 2; Outputs = 2
alpha = 0.5; Returns to Scale = VRS
─────────────────────────────────────────────────────────────
   efficiency     slackX1     slackX2     slackY1     slackY2
─────────────────────────────────────────────────────────────
1     0.68185  0.605935    4.26672e-8  3.91163e-8  4.67865
2     1.0      5.34772e-8  6.70059e-8  2.7034e-8   4.8232e-8
3     1.0      4.82929e-8  7.26916e-8  6.00225e-8  1.66806e-8
4     0.25     4.6491e-8   9.94558e-8  5.39305e-8  9.75587e-8
5     0.36     0.2         3.4         3.0         8.07052e-8
─────────────────────────────────────────────────────────────
```
"""
function deagdf(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    alpha::Float64 = 0.5, rts::Symbol = :CRS, slack::Bool = true,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix,Vector,Nothing} = nothing,
    names::Union{Vector{String},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::GeneralizedDFDEAModel

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    if Xref === nothing Xref = X end
    if Yref === nothing Yref = Y end

    nrefx, mref = size(Xref, 1), size(Xref, 2)
    nrefy, sref = size(Yref, 1), size(Yref, 2)

    if nx != ny
        throw(DimensionMismatch("number of rows in X and Y ($nx, $ny) are not equal"));
    end
    if nrefx != nrefy
        throw(DimensionMismatch("number of rows in Xref and Yref ($nrefx, $nrefy) are not equal"));
    end
    if m != mref
        throw(DimensionMismatch("number of columns in X and Xref ($m, $mref) are not equal"));
    end
    if s != sref
        throw(DimensionMismatch("number of columns in Y and Yref ($s, $sref) are not equal"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:NLP)
    end

    # Compute efficiency for each DMU
    n = nx
    nref = nrefx

    effi = zeros(n)
    lambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Create the optimization model
        deamodel = newdeamodel(optimizer)

        @variable(deamodel, eff, start = 1.0)
        @variable(deamodel, lambda[1:nref] >= 0)

        @NLobjective(deamodel, Min, eff)

        @NLconstraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= eff^(1-alpha) * x0[j])
        @NLconstraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= y0[j] / (eff^alpha) )

        # Add return to scale constraints
        if rts == :CRS
            # No contraint to add for constant returns to scale
        elseif rts == :VRS
            @constraint(deamodel, sum(lambda) == 1)
        else
            throw(ArgumentError("`rts` must be :CRS or :VRS"));
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        effi[i]  = JuMP.objective_value(deamodel)
        lambdaeff[i,:] = JuMP.value.(lambda)

        # Check termination status
        if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    # Get first-stage X and Y targets
    Xtarget = X .* effi .^(1-alpha)
    Ytarget = Y ./ ( effi .^alpha )

    # Compute slacks
    if slack == true
        # Use additive model with X and Y targets to get slacks
        slacksmodel = deaadd(Xtarget, Ytarget, :Ones, rts = rts, Xref = Xref, Yref = Yref, optimizer = optimizer)
        slackX = slacks(slacksmodel, :X)
        slackY = slacks(slacksmodel, :Y)

        # Get second-stage X and Y targets
        Xtarget = Xtarget - slackX
        Ytarget = Ytarget + slackY
    else
        if typeof(Xtarget) <: AbstractVector    Xtarget = Xtarget[:,:]  end
        if typeof(Ytarget) <: AbstractVector    Ytarget = Ytarget[:,:]  end

        slackX = Array{Float64}(undef, 0, 0)
        slackY = Array{Float64}(undef, 0, 0)
    end

    return GeneralizedDFDEAModel(n, m, s, alpha, rts, names, effi, slackX, slackY, lambdaeff, Xtarget, Ytarget)

end

function Base.show(io::IO, x::GeneralizedDFDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    eff = efficiency(x)
    dmunames = names(x)

    slackX = slacks(x, :X)
    slackY = slacks(x, :Y)
    hasslacks = ! isempty(slackX)

    if !compact
        print(io, "Generalized DF DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "alpha = ", x.alpha)
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        if hasslacks == true
            show(io, CoefTable(hcat(eff, slackX, slackY), ["efficiency"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
        else
            show(io, CoefTable(hcat(eff), ["efficiency"], dmunames))
        end
    end

end
