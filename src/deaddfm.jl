# This file contains functions for the Directional Multiplier DEA model
"""
    DirectionalMultiplierDEAModel
An data structure representing a directional distance function multiplier DEA model.
"""
struct DirectionalMultiplierDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64
    Gx::Symbol
    Gy::Symbol
    rts::Symbol
    dmunames::Union{Vector{AbstractString},Nothing}
    eff::Vector
    v::Matrix
    u::Matrix    
    omega::Vector
    Xtarget::Matrix
    Ytarget::Matrix
end

"""
    deaddfm(X, Y; Gx, Gy)
Compute data envelopment analysis directional distance function multiplier model for inputs
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
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `names`: a vector of strings with the names of the decision making units.
"""
function deaddfm(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    Gx::Union{Symbol, Matrix, Vector}, Gy::Union{Symbol, Matrix, Vector},
    rts::Symbol = :CRS,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix,Vector,Nothing} = nothing,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::DirectionalMultiplierDEAModel

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
    vi = zeros(n, m)
    ui = zeros(n, s)
    omegai = zeros(n)

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

            @variable(deamodel, v[1:m] >= 0)
            @variable(deamodel, u[1:s] >= 0)      
            @variable(deamodel, omega)  

            @objective(deamodel, Min, - sum(u[r] * y0[r] for r in 1:s) + sum(v[i] * x0[i] for i in 1:m) + omega)

            @constraint(deamodel, [j in 1:nref], sum(u[r] * Yref[j,r] for r in 1:s) - sum(v[i] * Xref[j,i] for i in 1:m) - omega <= 0)
            @constraint(deamodel, sum(u[r] * Gy0[r] for r in 1:s) + sum(v[i] * Gx0[i] for i in 1:m) == 1)

            # Add return to scale constraints
            if rts == :CRS
                @constraint(deamodel, omega == 0)
            elseif rts == :VRS            
                # No contraint to add for variable returns to scale
            else
                throw(ArgumentError("`rts` must be :CRS or :VRS"));
            end

            # Optimize and return results
            JuMP.optimize!(deamodel)

            effi[i]  = JuMP.objective_value(deamodel)
            vi[i, :] = JuMP.value.(v)
            ui[i, :] = JuMP.value.(u)
            omegai[i] = JuMP.value(omega)

            # Check termination status
            if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
                @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
            end
        else
            effi[i]  = 0.0
            vi[i,:] .= 0.0
            ui[i,:] .= 0.0
            omegai[i] = 0.0
        end

    end

    # Get first-stage X and Y targets
    Xtarget = X .- effi .* Gx
    Ytarget = Y .+ effi .* Gy

    if typeof(Xtarget) <: AbstractVector    Xtarget = Xtarget[:,:]  end
    if typeof(Ytarget) <: AbstractVector    Ytarget = Ytarget[:,:]  end

    return DirectionalMultiplierDEAModel(n, m, s, Gxsym, Gysym, rts, names, effi, vi, ui, omegai, Xtarget, Ytarget)

end

function Base.show(io::IO, x::DirectionalMultiplierDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    dmunames = names(x)

    eff = efficiency(x)
    v = multipliers(x, :X)
    u = multipliers(x, :Y)

    if !compact
        print(io, "Directional DF DEA Model (Multiplier form)\n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Returns to Scale = ", string(x.rts))
        print(io, "\n")
        print(io, "Gx = ", string(x.Gx), "; Gy = ", string(x.Gy))
        print(io, "\n")

        show(io, CoefTable(hcat(eff, v, u), ["efficiency"; ["v$i" for i in 1:m ]; ["u$i" for i in 1:s ]], dmunames))
    end
end
