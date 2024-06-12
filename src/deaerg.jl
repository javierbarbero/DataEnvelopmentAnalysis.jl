# This file contains functions for the Enhanced Russell Graph Slack Based Measure DEA model
"""
    EnhanceRussellGraphDEAModel
An data structure representing an Enhanced Russell Graph Slack Based Measure DEA model.
"""
struct EnhancedRussellGraphDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64
    rts::Symbol
    dmunames::Union{Vector{AbstractString},Nothing}
    eff::Vector
    beta::Vector
    tX::Matrix
    tY::Matrix
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    Xtarget::Matrix
    Ytarget::Matrix
end

"""
    deaerg(X, Y)
Compute data envelopment analysis Enhanced Russell Graph Slack Based Measure for inputs `X` and outputs `Y`.

# Optional Arguments
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `names`: a vector of strings with the names of the decision making units.
"""
function deaerg(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector}; 
    rts::Symbol = :CRS,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix,Vector,Nothing} = nothing,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::EnhancedRussellGraphDEAModel

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
        optimizer = DEAOptimizer(:LP)
    end

    # Compute efficiency for each DMU
    n = nx
    nref = nrefx

    effi = zeros(n)
    betai = zeros(n)
    tXi = zeros(n, m)
    tYi = zeros(n, s)
    lambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Create the optimization model
        deamodel = newdeamodel(optimizer)
        
        @variable(deamodel, beta >= 0)
        @variable(deamodel, tX[1:m] >= 0)
        @variable(deamodel, tY[1:s] >= 0)
        @variable(deamodel, mu[1:nref] >= 0)

        @objective(deamodel, Min, beta - 1 / m * sum(tX[t] / x0[t] for t in 1:m))

        @constraint(deamodel, beta + 1 / s * sum(tY[t] / y0[t] for t in 1:s) == 1)
        @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * mu[t] for t in 1:nref) == beta * x0[j] - tX[j])
        @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * mu[t] for t in 1:nref) == beta * y0[j] + tY[j])

        # Add return to scale constraints
        if rts == :CRS
            # No contraint to add for constant returns to scale
        elseif rts == :VRS
            @constraint(deamodel, sum(mu) == beta)
        elseif rts == :FDH
            @constraint(deamodel, sum(mu) == beta)

            # Add FDH constraints
            @variable(deamodel, lambda[1:nref] >= 0)
            set_binary.(lambda[1:nref])

            @constraint(deamodel, sum(lambda) == 1)
            @constraint(deamodel, [t in 1:nref], mu[t] == beta * lambda[t])
        else
            throw(ArgumentError("`rts` must be :CRS or :VRS"));
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        effi[i]  = JuMP.objective_value(deamodel)
        betai[i] = JuMP.value(beta)
        tXi[i,:] = JuMP.value.(tX)
        tYi[i,:] = JuMP.value.(tY)
        mui = JuMP.value.(mu)
        lambdaeff[i,:] = mui ./ betai[i]

        # Check termination status
        if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    slackX = tXi ./ betai
    slackY = tYi ./ betai

    # Get X and Y targets
    Xtarget = X - slackX
    Ytarget = Y + slackY

    return EnhancedRussellGraphDEAModel(n, m, s, rts, names, effi, betai, tXi, tYi, slackX, slackY, lambdaeff, Xtarget, Ytarget)

end

function Base.show(io::IO, x::EnhancedRussellGraphDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    dmunames = names(x)

    eff = efficiency(x)
    beta = x.beta
    slackX = slacks(x, :X)
    slackY = slacks(x, :Y)
    hasslacks = ! isempty(slackX)

    if !compact
        print(io, "Enhanced Russell Graph Slack Based Measure DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = Graph")
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        show(io, CoefTable(hcat(eff, beta, slackX, slackY), ["efficiency"; "beta"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
    end

end

function efficiency(model::EnhancedRussellGraphDEAModel, type::Symbol)::Vector

    if type == :beta return model.beta end

    throw(ArgumentError("$(typeof(model)) has no efficiency $(type)"));

end
