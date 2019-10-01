"""
    DirectionalDEAModel
An data structure representing a directional distance function DEA model.
"""
struct DirectionalDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64
    rts::Symbol
    eff::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
end

"""
    deaddf(X, Y, Gx, Gy)
Compute data envelopment analysis directional distance function model for inputs
`X` and outputs `Y`, using directions `Gx` and `Gy`.

# Optional Arguments
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `slack=true`: computes input and output slacks.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaddf(X, Y, ones(size(X)), ones(size(Y)))
Directional DF DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Returns to Scale = CRS
────────────────
      efficiency
────────────────
1   -3.43053e-16
2    3.21996
3    2.12169
4    0.0
5    6.73567
6    1.94595
7    0.0
8    3.63586
9    1.83784
10  10.2311
11   0.0
────────────────
```
"""
function deaddf(X::Matrix, Y::Matrix, Gx::Matrix, Gy::Matrix; rts::Symbol = :CRS, slack = true, Xref::Matrix = X, Yref::Matrix = Y)::DirectionalDEAModel
    # Check parameters
    nx, m = size(X)
    ny, s = size(Y)

    nGx, mGx = size(Gx)
    nGy, sGy = size(Gy)

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
    if size(Gx) != size(X)
        error("size of inputs should be equal to size of inputs direction")
    end
    if size(Gy) != size(Y)
        error("size of outputs should be equal to size of outputs direction")
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
        deamodel = Model(with_optimizer(GLPK.Optimizer))
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
            error("Invalid returns to scale $rts. Returns to scale should be :CRS or :VRS")
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        effi[i]  = JuMP.objective_value(deamodel)
        lambdaeff[i,:] = JuMP.value.(lambda)

    end

    # Compute slacks
    if slack == true

        # Get first-stage efficient X and Y
        Xeff = X .- effi .* Gx
        Yeff = Y .+ effi .* Gy

        # Use additive model with radial efficient X and Y to get slacks
        radialSlacks = deaadd(Xeff, Yeff, :Ones, rts = rts, Xref = Xref, Yref = Yref)
        slackX = slacks(radialSlacks, :X)
        slackY = slacks(radialSlacks, :Y)
    else
        slackX = Array{Float64}(undef, 0, 0)
        slackY = Array{Float64}(undef, 0, 0)
    end

    return DirectionalDEAModel(n, m, s, rts, effi, slackX, slackY, lambdaeff)

end

function deaddf(X::Vector, Y::Matrix, Gx::Vector, Gy::Matrix; rts::Symbol = :CRS, slack = true, Xref::Vector = X, Yref::Matrix = Y)::DirectionalDEAModel
    X = X[:,:]
    Xref = Xref[:,:]
    Gx = Gx[:,:]
    return deaddf(X, Y, Gx, Gy, rts = rts, slack = slack, Xref = Xref, Yref = Yref)
end

function deaddf(X::Matrix, Y::Vector, Gx::Matrix, Gy::Vector; rts::Symbol = :CRS, slack = true, Xref::Matrix = X, Yref::Vector = Y)::DirectionalDEAModel
    Y = Y[:,:]
    Yref = Yref[:,:]
    Gy = Gy[:,:]
    return deaddf(X, Y, Gx, Gy, rts = rts, slack = slack, Xref = Xref, Yref = Yref)
end

function deaddf(X::Vector, Y::Vector, Gx::Vector, Gy::Vector; rts::Symbol = :CRS, slack = true, Xref::Vector = X, Yref::Vector = Y)::DirectionalDEAModel
    X = X[:,:]
    Xref = Xref[:,:]
    Gx = Gx[:,:]
    Y = Y[:,:]
    Yref = Yref[:,:]
    Gy = Gy[:,:]
    return deaddf(X, Y, Gx, Gy, rts = rts, slack = slack, Xref = Xref, Yref = Yref)
end

function Base.show(io::IO, x::DirectionalDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    eff = efficiency(x)
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
        if hasslacks == true
            show(io, CoefTable(hcat(eff, slackX, slackY), ["efficiency"; ["slackX$i" for i in 1:m ]; ; ["slackY$i" for i in 1:s ]], ["$i" for i in 1:n]))
        else
            show(io, CoefTable(hcat(eff), ["efficiency"], ["$i" for i in 1:n]))
        end
    end

end
