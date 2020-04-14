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
    dmunames::Vector{String}
    eff::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
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

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> dea(X, Y)
Radial DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Orientation = Input; Returns to Scale = CRS
──────────────────────────────────────────────────
    efficiency       slackX1      slackX2  slackY1
──────────────────────────────────────────────────
1     1.0        0.0          0.0              0.0
2     0.62229   -4.41868e-15  0.0              0.0
3     0.819856   0.0          8.17926e-15      0.0
4     1.0       -8.03397e-16  0.0              0.0
5     0.310371   1.80764e-15  0.0              0.0
6     0.555556   4.44444      0.0              0.0
7     1.0        0.0          0.0              0.0
8     0.757669   1.60679e-15  0.0              0.0
9     0.820106   1.64021      0.0              0.0
10    0.490566   9.68683e-15  0.0              0.0
11    1.0        0.0          4.0              0.0
──────────────────────────────────────────────────
```
"""
function dea(X::Matrix, Y::Matrix; orient::Symbol = :Input, rts::Symbol = :CRS, slack = true, Xref::Matrix = X, Yref::Matrix = Y,
    disposX::Symbol = :Strong, disposY::Symbol = :Strong, 
    names::Vector{String} = Array{String}(undef, 0))::RadialDEAModel

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

    if disposX != :Strong && disposX != :Weak
        error("Invalued inputs disposability $disposX. Disposability should be :Strong or :Weak")
    end

    if disposY != :Strong && disposY != :Weak
        error("Invalid outputs disposability $disposY. Disposability should be :Strong or :Weak")
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
        deamodel = Model(GLPK.Optimizer)
        
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
            error("Invalid orientation $orient. Orientation should be :Input or :Output")
        end

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
        if termination_status(deamodel) != MOI.OPTIMAL
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    # Compute slacks
    if slack == true

        # Get first-stage efficient X and Y
        if orient == :Input
            Xeff = X .* effi
            Yeff = Y
        elseif orient == :Output
            Xeff = X
            Yeff = Y .* effi
        end

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

        radialSlacks = deaadd(Xeff, Yeff, rhoX = rhoX, rhoY = rhoY, rts = rts, Xref = Xref, Yref = Yref)
        slackX = slacks(radialSlacks, :X)
        slackY = slacks(radialSlacks, :Y)
    else
        slackX = Array{Float64}(undef, 0, 0)
        slackY = Array{Float64}(undef, 0, 0)
    end

    return RadialDEAModel(n, m, s, orient, rts, disposX, disposY, names, effi, slackX, slackY, lambdaeff)

end

function dea(X::Vector, Y::Matrix; orient::Symbol = :Input, rts::Symbol = :CRS, slack = true, Xref::Vector = X, Yref::Matrix = Y,
    disposX::Symbol = :Strong, disposY::Symbol = :Strong,
    names::Vector{String} = Array{String}(undef, 0))::RadialDEAModel

    X = X[:,:]
    Xref = Xref[:,:]
    return dea(X, Y, orient = orient, rts = rts, slack = slack, Xref = Xref, Yref = Yref, disposX = disposX, disposY = disposY, names = names)
end

function dea(X::Matrix, Y::Vector; orient::Symbol = :Input, rts::Symbol = :CRS, slack = true, Xref::Matrix = X, Yref::Vector = Y,
    disposX::Symbol = :Strong, disposY::Symbol = :Strong,
    names::Vector{String} = Array{String}(undef, 0))::RadialDEAModel

    Y = Y[:,:]
    Yref = Yref[:,:]
    return dea(X, Y, orient = orient, rts = rts, slack = slack, Xref = Xref, Yref = Yref, disposX = disposX, disposY = disposY, names = names)
end

function dea(X::Vector, Y::Vector; orient::Symbol = :Input, rts::Symbol = :CRS, slack = true, Xref::Vector = X, Yref::Vector = Y,
    disposX::Symbol = :Strong, disposY::Symbol = :Strong,
    names::Vector{String} = Array{String}(undef, 0))::RadialDEAModel

    X = X[:,:]
    Xref = Xref[:,:]
    Y = Y[:,:]
    Yref = Yref[:,:]
    return dea(X, Y, orient = orient, rts = rts, slack = slack, Xref = Xref, Yref = Yref, disposX = disposX, disposY = disposY, names = names)
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
        if disposX == :Weak println("Weak disposability of inputs") end
        if disposY == :Weak println("Weak disposability of outputs") end

        if hasslacks == true
            show(io, CoefTable(hcat(eff, slackX, slackY), ["efficiency"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
        else
            show(io, CoefTable(hcat(eff), ["efficiency"], dmunames))
        end
    end

end
