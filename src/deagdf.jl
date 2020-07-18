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
    dmunames::Vector{String}
    eff::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
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
2     1.0      5.19411e-8  6.85405e-8  2.54977e-8  4.97494e-8
3     1.0      4.5869e-8   7.51401e-8  5.75966e-8  1.92107e-8
4     0.25     4.71611e-8  9.87945e-8  5.46022e-8  9.69031e-8
5     0.36     0.2         3.4         3.0         8.07052e-8
─────────────────────────────────────────────────────────────
```
"""
function deagdf(X::Matrix, Y::Matrix; alpha::Float64 = 0.5, rts::Symbol = :CRS, slack = true, Xref::Matrix = X, Yref::Matrix = Y,
    names::Vector{String} = Array{String}(undef, 0))::GeneralizedDFDEAModel

    # Check parameters
    nx, m = size(X)
    ny, s = size(Y)

    nrefx, mref = size(Xref)
    nrefy, sref = size(Yref)

    if nx != ny
        error("number of observations is different in inputs and outputs")
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

    effi = zeros(n)
    lambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Create the optimization model
        deamodel = Model(Ipopt.Optimizer)
        set_silent(deamodel)
        set_optimizer_attribute(deamodel, "print_level", 0)

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
            error("Invalid returns to scale $rts. Returns to scale should be :CRS or :VRS")
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        effi[i]  = JuMP.objective_value(deamodel)
        lambdaeff[i,:] = JuMP.value.(lambda)

        # Check termination status
        if termination_status(deamodel) != MOI.LOCALLY_SOLVED
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    # Get first-stage efficient X and Y
    Xeff = X .* effi .^(1-alpha)
    Yeff = Y ./ ( effi .^alpha )

    # Compute slacks
    if slack == true

        slackX = zeros(n, m)
        slackY = zeros(n, s)
        lambdaeff = spzeros(n, nref)

        for i=1:n
            # Value of inputs and outputs to evaluate
            x0 = Xeff[i,:]
            y0 = Yeff[i,:]

            # Create the optimization model
            deamodel = Model(Ipopt.Optimizer)
            set_silent(deamodel)
            set_optimizer_attribute(deamodel, "print_level", 0)

            @variable(deamodel, sX[1:m] >= 0)
            @variable(deamodel, sY[1:s] >= 0)
            @variable(deamodel, lambda[1:nref] >= 0)

            @objective(deamodel, Max, sum(sX[j] for j in 1:m) + sum(sY[j] for j in 1:s) )

            @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) == x0[j] - sX[j])
            @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) == y0[j] + sY[j])

            # Add return to scale constraints
            if rts == :CRS
                # No contraint to add for constant returns to scale
            elseif rts == :VRS
                @constraint(deamodel, sum(lambda) == 1)
            end
            # No need for else statement as RTS parameter already checked in first stage

            # Optimize and return results
            JuMP.optimize!(deamodel)

            lambdaeff[i,:] = JuMP.value.(lambda)

            slackX[i,:] = JuMP.value.(sX)
            slackY[i,:] = JuMP.value.(sY)

            # Check termination status
            if termination_status(deamodel) != MOI.LOCALLY_SOLVED
                @warn ("DMU $i Slacks termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
            end

        end
    else
        slackX = Array{Float64}(undef, 0, 0)
        slackY = Array{Float64}(undef, 0, 0)
    end

    return GeneralizedDFDEAModel(n, m, s, alpha, rts, names, effi, slackX, slackY, lambdaeff)

end

function deagdf(X::Vector, Y::Matrix; alpha::Float64 = 0.5, rts::Symbol = :CRS, slack = true, Xref::Vector = X, Yref::Matrix = Y,
    names::Vector{String} = Array{String}(undef, 0))::GeneralizedDFDEAModel

    X = X[:,:]
    Xref = Xref[:,:]
    return deagdf(X, Y, alpha = alpha, rts = rts, slack = slack, Xref = Xref, Yref = Yref, names = names)
end

function deagdf(X::Matrix, Y::Vector; alpha::Float64 = 0.5, rts::Symbol = :CRS, slack = true, Xref::Matrix = X, Yref::Vector = Y,
    names::Vector{String} = Array{String}(undef, 0))::GeneralizedDFDEAModel

    Y = Y[:,:]
    Yref = Yref[:,:]
    return deagdf(X, Y, alpha = alpha, rts = rts, slack = slack, Xref = Xref, Yref = Yref, names = names)
end

function deagdf(X::Vector, Y::Vector; alpha::Float64 = 0.5, rts::Symbol = :CRS, slack = true, Xref::Vector = X, Yref::Vector = Y,
    names::Vector{String} = Array{String}(undef, 0))::GeneralizedDFDEAModel

    X = X[:,:]
    Xref = Xref[:,:]
    Y = Y[:,:]
    Yref = Yref[:,:]
    return deagdf(X, Y, alpha = alpha, rts = rts, slack = slack, Xref = Xref, Yref = Yref, names = names)
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
