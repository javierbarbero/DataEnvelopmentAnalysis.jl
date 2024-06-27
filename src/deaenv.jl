"""
    EnviromentalDEAModel
An data structure representing an enviromental DEA model.
"""
struct EnvironmentalDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64
    b::Int64
    Gx::Symbol
    Gy::Symbol
    Gb::Symbol
    rts::Symbol
    dmunames::Union{Vector{AbstractString},Nothing}
    eff::Vector
    slackX::Matrix
    slackY::Matrix
    slackB::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    Xtarget::Matrix
    Ytarget::Matrix
    Btarget::Matrix
end

"""
    deaenv(X, Y, B; Gx, Gy, Gb)
Compute data envelopment analysis environmental model for inputs
`X`, good outputs `Y`, and bad outputs `B`, using directions `Gx`, `Gy`, and `Gb``.

# Direction specification:

The directions `Gx`, `Gy`, and `Gb` can be one of the following symbols.
- `:Zeros`: use zeros.
- `:Ones`: use ones.
- `:Observed`: use observed values.
- `:Mean`: use column means.

Alternatively, a vector or matrix with the desired directions can be supplied.

# Optional Arguments
- `rts=:CRS`: chooses constant returns to scale.
- `slack=true`: computes input and output slacks.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of good outputs against which the units are evaluated.
- `Bref=B`: Identifies the reference set of bad outputs against which the units are evaluated.
- `names`: a vector of strings with the names of the decision making units.
"""
function deaenv(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector}, B::Union{Matrix,Vector};
    Gx::Union{Symbol, Matrix, Vector} = :Zeros, Gy::Union{Symbol, Matrix, Vector} = :Observed, Gb::Union{Symbol, Matrix, Vector} = :Observed,
    rts::Symbol = :CRS, slack::Bool = true,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix,Vector,Nothing} = nothing, Bref::Union{Matrix,Vector,Nothing} = nothing,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::EnvironmentalDEAModel

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)
    nb, b = size(B, 1), size(B, 2)

    if Xref === nothing Xref = X end
    if Yref === nothing Yref = Y end
    if Bref === nothing Bref = B end

    nrefx, mref = size(Xref, 1), size(Xref, 2)
    nrefy, sref = size(Yref, 1), size(Yref, 2)
    nrefb, bref = size(Bref, 1), size(Bref, 2)

    # Check data dimensions
    nx == ny || throw(DimensionMismatch("number of rows in X and Y ($nx, $ny) are not equal"));
    nrefx == nrefy || throw(DimensionMismatch("number of rows in Xref and Yref ($nrefx, $nrefy) are not equal"));
    m == mref || throw(DimensionMismatch("number of columns in X and Xref ($m, $mref) are not equal"));
    s == sref || throw(DimensionMismatch("number of columns in Y and Yref ($s, $sref) are not equal"));
    ny == nb || throw(DimensionMismatch("number of rows in Y and B ($ny, $nb) are not equal"));
    b == bref || throw(DimensionMismatch("number of columns in Y and Yref ($b, $bref) are not equal"));
    nrefy == nrefb || throw(DimensionMismatch("number of rows in Yref and Bref ($nrefy, $nrefb) are not equal"));

    # Build or get user directions
    if typeof(Gx) == Symbol
        Gxsym = Gx

        if Gx == :Zeros
            Gx = zeros(size(X))
        elseif Gx == :Ones
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

        if Gy == :Zeros
            Gy = zeros(size(Y))
        elseif Gy == :Ones
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

    if typeof(Gb) == Symbol
        Gbsym = Gb

        if Gb == :Zeros
            Gb = zeros(size(B))
        elseif Gb == :Ones
            Gb = ones(size(B))
        elseif Gb == :Observed
            Gb = B
        elseif Gb == :Mean
            Gb = repeat(mean(B, dims = 1), size(B, 1))
        else
            throw(ArgumentError("Invalid `Gb`"));
        end
    else
        Gbsym = :Custom
    end

    # Check directions dimensions
    nGx, mGx = size(Gx, 1), size(Gx, 2)
    nGy, sGy = size(Gy, 1), size(Gy, 2)
    nGb, bGb = size(Gy, 1), size(Gy, 2)

    (nGx == nx) & (mGx == m) || throw(DimensionMismatch("size of Gx and X ($(size(Gx)), $(size(X))) are not equal"));
    (nGy == ny) & (sGy == s) || throw(DimensionMismatch("size of Gy and Y ($(size(Gy)), $(size(Y))) are not equal"));
    (nGb == nb) & (bGb == b) || throw(DimensionMismatch("size of Gb and B ($(size(Gb)), $(size(B))) are not equal"));

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:LP)
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
        b0 = B[i,:]

        # Directions to use
        Gx0 = Gx[i,:]
        Gy0 = Gy[i,:]
        Gb0 = Gb[i,:]

        # Maximum of bad outputs
        maxB = maximum([maximum(B) maximum(Bref)])

        # Solve if any direction is different from zero
        if any(Gx0 .!= 0) | any(Gy0 .!= 0) | any(Gb .!= 0)
            # Create the optimization model
            deamodel = newdeamodel(optimizer)

            @variable(deamodel, eff)
            @variable(deamodel, lambda[1:nref] >= 0)

            @objective(deamodel, Max, eff)

            @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= x0[j] - eff * Gx0[j])
            @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= y0[j] + eff * Gy0[j])
            @constraint(deamodel, [j in 1:b], sum(Bref[t,j] * lambda[t] for t in 1:nref) <= b0[j] - eff * Gb0[j])
            @constraint(deamodel, [j in 1:b], maxB >= b0[j])

            # Add return to scale constraints
            if rts == :CRS
                # No contraint to add for constant returns to scale
            else
                throw(ArgumentError("`rts` must be :CRS"));
            end

            # Optimize and return results
            JuMP.optimize!(deamodel)

            effi[i]  = JuMP.objective_value(deamodel)
            lambdaeff[i,:] = JuMP.value.(lambda)

            # Check termination status
            if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
                @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
            end
        else
            effi[i]  = 0.0
            lambdaeff[i,:] .= 0.0
            lambdaeff[i,i] = 1.0
        end

    end

    # Get first-stage X and Y targets
    Xtarget = X .- effi .* Gx
    Ytarget = Y .+ effi .* Gy
    Btarget = B .- effi .* Gb

    # Compute slacks
    if slack == true
        # Use additive model with X and Y targets to get slacks
        # Bad outputs are included as additional inputs.
        slacksmodel = deaadd([Xtarget Btarget], Ytarget, :Ones, rts = rts, Xref = [Xref Bref], Yref = Yref, optimizer = optimizer)
        slackX = slacks(slacksmodel, :X)[:,1:m]
        slackB = slacks(slacksmodel, :X)[:,(m+1):(m+b)]
        slackY = slacks(slacksmodel, :Y)

        # Get second-stage X and Y targets
        Xtarget = Xtarget - slackX
        Ytarget = Ytarget + slackY
        Btarget = Btarget - slackB
    else
        if typeof(Xtarget) <: AbstractVector    Xtarget = Xtarget[:,:]  end
        if typeof(Ytarget) <: AbstractVector    Ytarget = Ytarget[:,:]  end
        if typeof(Btarget) <: AbstractVector    Btarget = Btarget[:,:]  end

        slackX = Array{Float64}(undef, 0, 0)
        slackY = Array{Float64}(undef, 0, 0)
        slackB = Array{Float64}(undef, 0, 0)
    end

    return EnvironmentalDEAModel(n, m, s, b, Gxsym, Gysym, Gbsym, rts, names, effi, slackX, slackY, slackB, lambdaeff, Xtarget, Ytarget, Btarget)

end

isenvironmental(::EnvironmentalDEAModel) = true;

function Base.show(io::IO, x::EnvironmentalDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    b = nbadoutputs(x)
    eff = efficiency(x)
    dmunames = names(x)

    slackX = slacks(x, :X)
    slackY = slacks(x, :Y)
    slackB = slacks(x, :B)
    hasslacks = ! isempty(slackX)

    if !compact
        print(io, "Environmental DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "; Bad Outputs = ", b)
        print(io, "\n")
        print(io, "Returns to Scale = ", string(x.rts))
        print(io, "\n")
        print(io, "Gx = ", string(x.Gx), "; Gy = ", string(x.Gy), " ; Gb = ", string(x.Gb))
        print(io, "\n")

        if hasslacks == true
            show(io, CoefTable(hcat(eff, slackX, slackY, slackB), ["efficiency"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]; ["slackB$i" for i in 1:b ]], dmunames))
        else
            show(io, CoefTable(hcat(eff), ["efficiency"], dmunames))
        end
    end

end

