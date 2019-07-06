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
    eff::Vector
    lambda::SparseMatrixCSC{Float64, Int64}
end

"""
    dea(X, Y)
Computes data envelopment analysis radial model for inputs `X` and outputs `Y`.

### Optional Arguments
- `orient=:Input`: chosse between input oriented radial model `:Input` or
output oriented radial model `:Output`.
- `rts=:CRS`: chosse between constant returns to scale `:CRS` or variable
returns to scale `:VRS`.
- `Xref=X`: reference set of inputs to which evaluate the units.
- `Yref=Y`: reference set of outputs to which evaluate the units.

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];
julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];
julia> dea(X, Y)
Radial DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Orientation = Input; Returns to Scale = CRS
──────────────
    efficiency
──────────────
1     1.0
2     0.62229
3     0.819856
4     1.0
5     0.310371
6     0.555556
7     1.0
8     0.757669
9     0.820106
10    0.490566
11    1.0
──────────────
```
"""
function dea(X::Matrix, Y::Matrix; orient::Symbol = :Input, rts::Symbol = :CRS, Xref::Matrix = X, Yref::Matrix = Y)::RadialDEAModel
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

    efficiency = zeros(n)
    lambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Create the optimization model
        deamodel = Model(with_optimizer(GLPK.Optimizer))
        @variable(deamodel, eff)
        @variable(deamodel, lambda[1:nref] >= 0)

        if orient == :Input
            # Input orientation
            @objective(deamodel, Min, eff)

            @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= eff * x0[j])
            @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= y0[j])
        elseif orient == :Output
            # Output orientation
            @objective(deamodel, Max, eff)

            @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= x0[j])
            @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= eff * y0[j])
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

        efficiency[i]  = JuMP.objective_value(deamodel)
        lambdaeff[i,:] = JuMP.value.(lambda)

    end

    return RadialDEAModel(n, m, s, orient, rts, efficiency, lambdaeff)

end

function dea(X::Vector, Y::Matrix; orient::Symbol = :Input, rts::Symbol = :CRS, Xref::Vector = X, Yref::Matrix = Y)::RadialDEAModel
    X = X[:,:]
    Xref = X[:,:]
    return dea(X, Y, orient = orient, rts = rts, Xref = Xref, Yref = Yref)
end

function dea(X::Matrix, Y::Vector; orient::Symbol = :Input, rts::Symbol = :CRS, Xref::Matrix = X, Yref::Vector = Y)::RadialDEAModel
    Y = Y[:,:]
    Yref = Y[:,:]
    return dea(X, Y, orient = orient, rts = rts, Xref = Xref, Yref = Yref)
end

function dea(X::Vector, Y::Vector; orient::Symbol = :Input, rts::Symbol = :CRS, Xref::Vector = X, Yref::Vector = Y)::RadialDEAModel
    X = X[:,:]
    Xref = X[:,:]
    Y = Y[:,:]
    Yref = Y[:,:]
    return dea(X, Y, orient = orient, rts = rts, Xref = Xref, Yref = Yref)
end

function Base.show(io::IO, x::RadialDEAModel)
    compact = get(io, :compact, false)

    eff = efficiency(x)
    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)

    if !compact
        print(io, "Radial DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = ", string(x.orient))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        show(io, CoefTable(hcat(eff), ["efficiency"], ["$i" for i in 1:n]))
    else

    end
end
