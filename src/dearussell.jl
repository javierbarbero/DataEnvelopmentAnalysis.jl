# This file contains functions for the Russell DEA model
"""
    RussellDEAModel
An data structure representing a Russell DEA model.
"""
struct RussellDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64
    orient::Symbol
    rts::Symbol
    dmunames::Union{Vector{AbstractString},Nothing}
    eff::Vector
    thetaX::Matrix
    thetaY::Matrix
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    Xtarget::Matrix
    Ytarget::Matrix
end

"""
    dearussell(X, Y)
Compute the Russell model using data envelopment analysis for inputs X and outputs Y.

# Optional Arguments
- `orient=:Input`: chooses the Russell input mode. For the Russell output model choose `:Output`. For the Russell graph model choose `:Graph`.
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `slack=true`: computes input and output slacks.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `names`: a vector of strings with the names of the decision making units.
"""
function dearussell(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    orient::Symbol = :Input, rts::Symbol = :CRS, slack::Bool = true,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix,Vector,Nothing} = nothing,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::RussellDEAModel

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
        if orient == :Input || orient == :Output
            optimizer = DEAOptimizer(:LP)
        else
            optimizer = DEAOptimizer(:NLP)

            rts != :FDH || throw(ArgumentError("Free Disposal Hull, rts = :FDH, not implemented for orientation :Graph"))
        end
    end

    # Compute efficiency for each DMU
    n = nx
    nref = nrefx

    effi = zeros(n)
    thetaXi = zeros(n, m)
    thetaYi = zeros(n, s)
    lambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Create the optimization model
        if orient == :Input
            # Input orientation
            deamodel = newdeamodel(optimizer)
            
            mno0 = sum(x0 .!= 0)

            @variable(deamodel, lambda[1:nref] >= 0)
            @variable(deamodel, theta[1:m] <= 1)
            @objective(deamodel, Min, 1 / mno0 * sum(theta[t] for t in 1:m if x0[t] != 0 ))

            @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= theta[j] * x0[j])
            @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= y0[j])
        elseif orient == :Output
            # Output orientation
            deamodel = newdeamodel(optimizer)
            
            sno0 = sum(y0 .!= 0)

            @variable(deamodel, lambda[1:nref] >= 0)
            @variable(deamodel, theta[1:s] >= 1)
            @objective(deamodel, Max, 1 / sno0 * sum(theta[t] for t in 1:s if y0[t] != 0 ))

            @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= x0[j])
            @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= theta[j] * y0[j])
        elseif orient == :Graph
            # Graph orientation
            deamodel = newdeamodel(optimizer)
            set_silent(deamodel)            

            mno0 = sum(x0 .!= 0)
            sno0 = sum(y0 .!= 0)

            @variable(deamodel, lambda[1:nref] >= 0)
            @variable(deamodel, theta[1:m] <= 1)
            @variable(deamodel, phi[1:s] >= 1)

            @NLobjective(deamodel, Min, 1 / (mno0 + sno0) * (sum(theta[t] for t in 1:m if x0[t] != 0 ) + sum(1/phi[t] for t in 1:s if y0[t] != 0) ))

            @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) == theta[j] * x0[j])
            @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) == phi[j] * y0[j])
        else
            throw(ArgumentError("`orient` must be :Input, :Output or :Graph"));
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
            throw(ArgumentError("`rts` must be :CRS, :VRS or :FDH"));
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        effi[i]  = JuMP.objective_value(deamodel)
        lambdaeff[i,:] = JuMP.value.(lambda)

        if orient == :Input
            thetaXi[i,:] = JuMP.value.(theta)
        elseif orient == :Output
            thetaYi[i,:] = JuMP.value.(theta)
        elseif orient == :Graph
            thetaXi[i,:] = JuMP.value.(theta)
            thetaYi[i,:] = JuMP.value.(phi)
        end

        # Check termination status
        if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    # Get first-stage X and Y targets
    if orient == :Input
        Xtarget = X .* thetaXi
        Ytarget = Y
    elseif orient == :Output
        Xtarget = X
        Ytarget = Y .* thetaYi
    elseif orient == :Graph
        Xtarget = X .* thetaXi
        Ytarget = Y .* thetaYi
    end

    # Compute slacks
    if (slack == true) && (orient != :Graph)

        # Use additive model with Russell efficient X and Y to get slacks
        russellSlacks = deaadd(Xtarget, Ytarget, :Ones, rts = rts, Xref = Xref, Yref = Yref, optimizer = optimizer)
        slackX = slacks(russellSlacks, :X)
        slackY = slacks(russellSlacks, :Y)

        # Get second-stage X and Y targets
        Xtarget = Xtarget - slackX
        Ytarget = Ytarget + slackY
    else
        if typeof(Xtarget) <: AbstractVector    Xtarget = Xtarget[:,:]  end
        if typeof(Ytarget) <: AbstractVector    Ytarget = Ytarget[:,:]  end

        slackX = Array{Float64}(undef, 0, 0)
        slackY = Array{Float64}(undef, 0, 0)
    end

    if orient == :Input
        thetaYi = Array{Float64}(undef, 0, 0)
    elseif orient == :Output
        thetaXi = Array{Float64}(undef, 0, 0)
    end

    return RussellDEAModel(n, m, s, orient, rts, names, effi, thetaXi, thetaYi, slackX, slackY, lambdaeff, Xtarget, Ytarget)

end

function Base.show(io::IO, x::RussellDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    dmunames = names(x)

    eff = efficiency(x)
    slackX = slacks(x, :X)
    slackY = slacks(x, :Y)
    hasslacks = ! isempty(slackX)

    if !compact
        print(io, "Russell DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = ", string(x.orient))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        if x.orient == :Input
            thetaX = x.thetaX
            if hasslacks == true
                show(io, MIME"text/plain"(), CoefTable(hcat(eff, thetaX, slackY), ["efficiency"; ["effX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
            else
                show(io, MIME"text/plain"(), CoefTable(hcat(eff, thetaX), ["efficiency"; ["effX$i" for i in 1:m ] ], dmunames))
            end
        elseif x.orient == :Output
            thetaY = x.thetaY
            if hasslacks == true
                show(io, MIME"text/plain"(), CoefTable(hcat(eff, thetaY, slackX), ["efficiency"; ["effY$i" for i in 1:s ]; ["slackX$i" for i in 1:m ]], dmunames))
            else
                show(io, MIME"text/plain"(), CoefTable(hcat(eff, thetaY), ["efficiency"; ["effY$i" for i in 1:s ] ], dmunames))
            end
        elseif x.orient == :Graph
            thetaX = x.thetaX
            thetaY = x.thetaY
            show(io, MIME"text/plain"(), CoefTable(hcat(eff, thetaX, thetaY), ["efficiency"; ["effX$i" for i in 1:m ]; ["effY$i" for i in 1:s ] ], dmunames))
        end
    end

end

function efficiency(model::RussellDEAModel, type::Symbol)::Matrix

    if (type == :X && (model.orient == :Input  || model.orient == :Graph)) return model.thetaX end
    if (type == :Y && (model.orient == :Output || model.orient == :Graph)) return model.thetaY end
        
    throw(ArgumentError("$(typeof(model)) with orienation $(model.orient) has no efficiency $(type)"));

end
