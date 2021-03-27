"""
    DirectionalDEAModel
An data structure representing a directional distance function DEA model.
"""
struct DirectionalDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64
    Gx::Symbol
    Gy::Symbol
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
    deaddf(X, Y; Gx, Gy)
Compute data envelopment analysis directional distance function model for inputs
`X` and outputs `Y`, using directions `Gx` and `Gy`.

# Direction specification:

The directions `Gx` and `Gy` can be one of the following symbols.
- `:Zeros`: use zeros.
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

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaddf(X, Y, Gx = :Ones, Gy = :Ones)
Directional DF DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Returns to Scale = CRS
Gx = Ones; Gy = Ones
─────────────────────────────────────────────────────
      efficiency       slackX1       slackX2  slackY1
─────────────────────────────────────────────────────
1   -3.43053e-16   0.0           0.0              0.0
2    3.21996      -3.21359e-15   0.0              0.0
3    2.12169       0.0          -4.80367e-15      0.0
4    0.0          -8.03397e-16   0.0              0.0
5    6.73567      -2.41019e-15   0.0              0.0
6    1.94595      10.9189        0.0              0.0
7    0.0           0.0           0.0              0.0
8    3.63586       6.42718e-15   0.0              0.0
9    1.83784       4.75676       0.0              0.0
10  10.2311        6.12173e-15   0.0              0.0
11   0.0           0.0           4.0              0.0
─────────────────────────────────────────────────────
```
"""
function deaddf(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    Gx::Union{Symbol, Matrix, Vector}, Gy::Union{Symbol, Matrix, Vector},
    rts::Symbol = :CRS, slack::Bool = true,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix,Vector,Nothing} = nothing,
    names::Union{Vector{String},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::DirectionalDEAModel

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
        optimizer = DEAOptimizer(GLPK.Optimizer)
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

        # Directions to use
        Gx0 = Gx[i,:]
        Gy0 = Gy[i,:]

        # Create the optimization model
        deamodel = newdeamodel(optimizer)

        @variable(deamodel, eff)
        @variable(deamodel, lambda[1:nref] >= 0)

        @objective(deamodel, Max, eff)

        @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= x0[j] - eff * Gx0[j])
        @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= y0[j] + eff * Gy0[j])

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
    Xtarget = X .- effi .* Gx
    Ytarget = Y .+ effi .* Gy

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

    return DirectionalDEAModel(n, m, s, Gxsym, Gysym, rts, names, effi, slackX, slackY, lambdaeff, Xtarget, Ytarget)

end

function Base.show(io::IO, x::DirectionalDEAModel)
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
        print(io, "Directional DF DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Returns to Scale = ", string(x.rts))
        print(io, "\n")
        print(io, "Gx = ", string(x.Gx), "; Gy = ", string(x.Gy))
        print(io, "\n")

        if hasslacks == true
            show(io, CoefTable(hcat(eff, slackX, slackY), ["efficiency"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
        else
            show(io, CoefTable(hcat(eff), ["efficiency"], dmunames))
        end
    end

end
