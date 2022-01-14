# This file contains functions for the Radial Multiplier DEA model
"""
    AbstractRadialMultiplierDEAModel
An abstract type representing a radial DEA model.
"""
abstract type AbstractRadialMultiplierDEAModel <: AbstractTechnicalDEAModel end

"""
    RadialMultiplierDEAModel
An data structure representing a radial multiplier DEA model.
"""
struct RadialMultiplierDEAModel <: AbstractRadialMultiplierDEAModel
    n::Int64
    m::Int64
    s::Int64
    orient::Symbol
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
    deam(X, Y)
Compute the radial multiplier model using data envelopment analysis for inputs X and outputs Y.

# Optional Arguments
- `orient=:Input`: chooses the radially oriented input mode. For the radially oriented output model choose `:Output`.
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `names`: a vector of strings with the names of the decision making units.
"""
function deam(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    orient::Symbol = :Input, rts::Symbol = :CRS,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix, Vector,Nothing} = nothing,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing, progress::Bool = false)::RadialMultiplierDEAModel

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
    vi = zeros(n, m)
    ui = zeros(n, s)
    omegai = zeros(n)

    p = progressbarmeter(n, desc = "Computing Radial DEA Model...", progress = progress)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Create the optimization model
        deamodel = newdeamodel(optimizer)

        @variable(deamodel, v[1:m] >= 0)
        @variable(deamodel, u[1:s] >= 0)      
        @variable(deamodel, omega)  

        if orient == :Input
            @objective(deamodel, Max, sum(u[r] * y0[r] for r in 1:s) - omega)

            @constraint(deamodel, sum(v[i] * x0[i] for i in 1:m) == 1)
            @constraint(deamodel, [j in 1:nref], sum(u[r] * Yref[j,r] for r in 1:s) - sum(v[i] * Xref[j,i] for i in 1:m) - omega <= 0)

        elseif orient == :Output
            @objective(deamodel, Min, sum(v[i] * x0[i] for i in 1:m) + omega)

            @constraint(deamodel, sum(u[r] * y0[r] for r in 1:s) == 1)
            @constraint(deamodel, [j in 1:nref], sum(v[i] * Xref[j,i] for i in 1:m) - sum(u[r] * Yref[j,r] for r in 1:s) + omega >= 0)
        else
            throw(ArgumentError("`orient` must be :Input or :Output"));
        end

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

    if typeof(Xtarget) <: AbstractVector    Xtarget = Xtarget[:,:]  end
    if typeof(Ytarget) <: AbstractVector    Ytarget = Ytarget[:,:]  end

    return RadialMultiplierDEAModel(n, m, s, orient, rts, names, effi, vi, ui, omegai, Xtarget, Ytarget)
end

function Base.show(io::IO, x::RadialMultiplierDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    dmunames = names(x)

    eff = efficiency(x)
    v = multipliers(x, :X)
    u = multipliers(x, :Y)

    if !compact
        print(io, "Radial DEA Model (Multiplier form)\n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = ", string(x.orient))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")

        show(io, CoefTable(hcat(eff, v, u), ["efficiency"; ["v$i" for i in 1:m ]; ["u$i" for i in 1:s ]], dmunames))
    end
end
