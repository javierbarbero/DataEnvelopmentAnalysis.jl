# This file contains functions for the Modified Directional Distance Function DEA model
"""
    ModifiedDDFDEAModel
An data structure representing a modified directional distance function DEA model.
"""
struct ModifiedDDFDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64
    Gx::Symbol
    Gy::Symbol
    rts::Symbol
    dmunames::Union{Vector{AbstractString},Nothing}
    eff::Vector
    betax::Vector
    betay::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    Xtarget::Matrix
    Ytarget::Matrix
end

"""
    deamddf(X, Y; Gx, Gy)
Compute data envelopment analysis modified directional distance function model for inputs
`X` and outputs `Y`, using directions `Gx` and `Gy`.

# Direction specification:

The directions `Gx` and `Gy` can be one of the following symbols.
- `:Ones`: use ones.
- `:Observed`: use observed values.
- `:Mean`: use column means.

Alternatively, a vector or matrix with the desired directions can be supplied.

# Optional Arguments
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `slack=true`: computes input and output slacks.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `names`: a vector of strings with the names of the decision making units.
"""
function deamddf(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    Gx::Union{Symbol, Matrix, Vector}, Gy::Union{Symbol, Matrix, Vector},
    rts::Symbol = :CRS, slack::Bool = true,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix,Vector,Nothing} = nothing,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::ModifiedDDFDEAModel

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

    # Build or get user directions
    if typeof(Gx) == Symbol
        Gxsym = Gx

        if Gx == :Ones
            Gx = ones(size(X))
        elseif Gx == :Observed
            Gx = X
        elseif Gx == :Mean
            Gx = repeat(mean(X, dims = 1), size(X, 1))
        else
            throw(ArgumentError("Invalid `Gx`"));
        end

    else
        Gxsym = :Custom
    end

    if typeof(Gy) == Symbol
        Gysym = Gy

        if Gy == :Ones
            Gy = ones(size(Y))
        elseif Gy == :Observed
            Gy = Y
        elseif Gy == :Mean
            Gy = repeat(mean(Y, dims = 1), size(Y, 1))
        else
            throw(ArgumentError("Invalid `Gy`"));
        end

    else
        Gysym = :Custom
    end

    nGx, mGx = size(Gx, 1), size(Gx, 2)
    nGy, sGy = size(Gy, 1), size(Gy, 2)

    if (size(Gx, 1) != size(X, 1)) | (size(Gx, 2) != size(X, 2))
        throw(DimensionMismatch("size of Gx and X ($(size(Gx)), $(size(X))) are not equal"));
    end
    if (size(Gy, 1) != size(Y, 1)) | (size(Gy, 2) != size(Y, 2))
        throw(DimensionMismatch("size of Gy and Y ($(size(Gy)), $(size(Y))) are not equal"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:LP)
    end

    # Compute efficiency for each DMU
    n = nx
    nref = nrefx

    effi = zeros(n)
    betaxi = zeros(n)
    betayi = zeros(n)
    lambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Directions to use
        Gx0 = Gx[i,:]
        Gy0 = Gy[i,:]

        # Solve if any direction is different from zero
        if any(Gx0 .!= 0) | any(Gy0 .!= 0)
            # Create the optimization model
            deamodel = newdeamodel(optimizer)
            set_silent(deamodel)

            @variable(deamodel, betax >= 0)
            @variable(deamodel, betay >= 0)
            @variable(deamodel, lambda[1:nref] >= 0)

            @objective(deamodel, Max, betax + betay)

            @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= x0[j] - betax * Gx0[j])
            @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= y0[j] + betay * Gy0[j])

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
            betaxi[i] = JuMP.value(betax)
            betayi[i] = JuMP.value(betay)
            lambdaeff[i,:] = JuMP.value.(lambda)

            # Check termination status
            if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
                @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
            end
        else
            effi[i]  = 0.0
            betaxi[i] = 0.0
            betayi[i] = 0.0
            lambdaeff[i,:] .= 0.0
            lambdaeff[i,i] = 1.0
        end

    end

    # Get first-stage X and Y targets
    Xtarget = X .- betaxi .* Gx
    Ytarget = Y .+ betayi .* Gy

    # Compute slacks
    if slack == true
        # Compute slacks if more than one input or more than 1 output
        if (m > 1) || (s > 1)
            # Use additive model with X and Y targets to get slacks
            slacksmodel = deaadd(Xtarget, Ytarget, :Ones, rts = rts, Xref = Xref, Yref = Yref, optimizer = optimizer)
            slackX = slacks(slacksmodel, :X)
            slackY = slacks(slacksmodel, :Y)

            # Get second-stage X and Y targets
            Xtarget = Xtarget - slackX
            Ytarget = Ytarget + slackY
        else
            slackX = zeros(n, m)
            slackY = zeros(n, s)

            if typeof(Xtarget) <: AbstractVector    Xtarget = Xtarget[:,:]  end
            if typeof(Ytarget) <: AbstractVector    Ytarget = Ytarget[:,:]  end
        end
    else
        if typeof(Xtarget) <: AbstractVector    Xtarget = Xtarget[:,:]  end
        if typeof(Ytarget) <: AbstractVector    Ytarget = Ytarget[:,:]  end

        slackX = Array{Float64}(undef, 0, 0)
        slackY = Array{Float64}(undef, 0, 0)
    end

    return ModifiedDDFDEAModel(n, m, s, Gxsym, Gysym, rts, names, effi, betaxi, betayi, slackX, slackY, lambdaeff, Xtarget, Ytarget)

end

function Base.show(io::IO, x::ModifiedDDFDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    eff = efficiency(x)
    dmunames = names(x)

    betax = x.betax
    betay = x.betay
    slackX = slacks(x, :X)
    slackY = slacks(x, :Y)
    hasslacks = ! isempty(slackX)

    if !compact
        print(io, "Modified DDF DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Returns to Scale = ", string(x.rts))
        print(io, "\n")
        print(io, "Gx = ", string(x.Gx), "; Gy = ", string(x.Gy))
        print(io, "\n")

        if hasslacks == true
            show(io, MIME"text/plain"(), CoefTable(hcat(eff, betax, betay, slackX, slackY), ["efficiency"; "βx"; "βy"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
        else
            show(io, MIME"text/plain"(), CoefTable(hcat(eff, betax, betay), ["efficiency"; "βx"; "βy"], dmunames))
        end
    end

end

function efficiency(model::ModifiedDDFDEAModel, type::Symbol)::Vector

    if type == :X return model.betax end
    if type == :Y return model.betay end

    throw(ArgumentError("$(typeof(model)) has no efficiency $(type)"));

end
