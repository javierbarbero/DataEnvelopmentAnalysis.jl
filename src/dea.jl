# This file contains functions for the Radial DEA model
"""
    AbstractRadialDEAModel
An abstract type representing a radial DEA model.
"""
abstract type AbstractRadialDEAModel <: AbstractTechnicalDEAModel end

"""
    RadialDEAModel
An data structure representing a radial DEA model.
"""
struct RadialDEAModel <: AbstractRadialDEAModel
    n::Int64
    m::Int64
    s::Int64
    orient::Symbol
    rts::Symbol
    disposX::Symbol
    disposY::Symbol
    dmunames::Union{Vector{AbstractString},Nothing}
    eff::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    Xtarget::Matrix
    Ytarget::Matrix
end

"""
    dea(X, Y)
Compute the radial model using data envelopment analysis for inputs X and outputs Y.

# Optional Arguments
- `orient=:Input`: chooses the radially oriented input mode. For the radially oriented output model choose `:Output`.
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `slack=true`: computes input and output slacks.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `disposX=:Strong`: chooses strong disposability of inputs. For weak disposability choose `:Weak`.
- `disposY=:Strong`: chooses strong disposability of outputs. For weak disposability choose `:Weak`.
- `names`: a vector of strings with the names of the decision making units.
"""
function dea(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    orient::Symbol = :Input, rts::Symbol = :CRS, slack::Bool = true,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix, Vector,Nothing} = nothing,
    disposX::Symbol = :Strong, disposY::Symbol = :Strong,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing, progress::Bool = false)::RadialDEAModel

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

    if disposX != :Strong && disposX != :Weak
        throw(ArgumentError("`disposX` must be :Strong or :Weak"));
    end

    if disposY != :Strong && disposY != :Weak
        throw(ArgumentError("`disposY` must be :Strong or :Weak"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:LP)
    end

    # Compute efficiency for each DMU
    n = nx
    nref = nrefx

    effi = zeros(n)
    lambdaeff = spzeros(n, nref)

    p = progressbarmeter(n, desc = "Computing Radial DEA Model...", progress = progress)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Create the optimization model
        deamodel = newdeamodel(optimizer)

        @variable(deamodel, eff)
        @variable(deamodel, lambda[1:nref] >= 0)

        if orient == :Input
            # Input orientation
            @objective(deamodel, Min, eff)

            # Inequality or equality restrictions based on disposability
            if disposX == :Strong
                @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= eff * x0[j])
            elseif disposX == :Weak
                @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) == eff * x0[j])
            end

            if disposY == :Strong
                @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= y0[j])
            elseif disposY == :Weak
                @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) == y0[j])
            end

        elseif orient == :Output
            # Output orientation
            @objective(deamodel, Max, eff)

            # Inequality or equality restrictions based on disposability
            if disposX == :Strong
                @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= x0[j])
            elseif disposX == :Weak
                @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) == x0[j])
            end

            if disposY == :Strong
                @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= eff * y0[j])
            elseif disposY == :Weak
                @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) == eff * y0[j])
            end

        else
            throw(ArgumentError("`orient` must be :Input or :Output"));
        end

        # Add return to scale constraints
        if rts == :CRS
            # No contraint to add for constant returns to scale
        elseif rts == :VRS
            @constraint(deamodel, sum(lambda) == 1)
        elseif rts == :FDH
            @constraint(deamodel, sum(lambda) == 1)
            set_binary.(lambda[1:nref])
        else
            throw(ArgumentError("`rts` must be :CRS, :VRS, or :FDH"));
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        effi[i]  = JuMP.objective_value(deamodel)
        lambdaeff[i,:] = JuMP.value.(lambda)

        # Check termination status
        if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

        next!(p)
    end

    # Get first-stage X and Y targets
    if orient == :Input
        Xtarget = X .* effi
        Ytarget = Y
    elseif orient == :Output
        Xtarget = X
        Ytarget = Y .* effi
    end

    # Compute slacks
    if slack == true

        # Use additive model with radial efficient X and Y to get slacks
        if disposX == :Strong
            rhoX = ones(size(X))
        elseif disposX == :Weak
            rhoX = zeros(size(X))
        end

        if disposY == :Strong
            rhoY = ones(size(Y))
        elseif disposY == :Weak
            rhoY = zeros(size(Y))
        end

        slacksmodel = deaadd(Xtarget, Ytarget, rhoX = rhoX, rhoY = rhoY, rts = rts, Xref = Xref, Yref = Yref, optimizer = optimizer)
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

    return RadialDEAModel(n, m, s, orient, rts, disposX, disposY, names, effi, slackX, slackY, lambdaeff, Xtarget, Ytarget)

end

function Base.show(io::IO, x::RadialDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    disposX = x.disposX
    disposY = x.disposY
    dmunames = names(x)

    eff = efficiency(x)
    slackX = slacks(x, :X)
    slackY = slacks(x, :Y)
    hasslacks = ! isempty(slackX)

    if !compact
        print(io, "Radial DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = ", string(x.orient))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        if disposX == :Weak print(io, "Weak disposability of inputs \n") end
        if disposY == :Weak print(io, "Weak disposability of outputs \n") end

        if hasslacks == true
            show(io, CoefTable(hcat(eff, slackX, slackY), ["efficiency"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
        else
            show(io, CoefTable(hcat(eff), ["efficiency"], dmunames))
        end
    end

end
