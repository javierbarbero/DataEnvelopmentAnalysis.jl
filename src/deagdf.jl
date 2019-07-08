# This file contains functions for the Generalized Distance Function DEA model
"""
    GeneralizedDFDEAModel
An data structure representing a generalized distance function DEA model.
"""
struct GeneralizedDFDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64
    alpha::Float64
    rts::Symbol
    eff::Vector
    lambda::SparseMatrixCSC{Float64, Int64}
end

"""
    deagdf(X, Y, alpha)
Compute generalized distance function data envelopment analysis model for
inputs `X`, outputs `Y`, and `alpha`.

# Optional Arguments
- `rts=:CRS`: chosse between constant returns to scale `:CRS` or variable
returns to scale `:VRS`.
- `Xref=X`: reference set of inputs to which evaluate the units.
- `Yref=Y`: reference set of outputs to which evaluate the units.

# Examples
```jldoctest
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9];

julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6];

julia> deagdf(X, Y, 0.5, rts = :VRS)
Generalized DF DEA Model
DMUs = 5; Inputs = 2; Outputs = 2
alpha = 0.5; Returns to Scale = VRS
─────────────
   efficiency
─────────────
1     0.68185
2     1.0
3     1.0
4     0.25
5     0.36
─────────────
```
"""
function deagdf(X::Matrix, Y::Matrix, alpha::Float64; rts::Symbol = :CRS, Xref::Matrix = X, Yref::Matrix = Y)::GeneralizedDFDEAModel
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
        deamodel = Model(with_optimizer(Ipopt.Optimizer, print_level= 0 ))
        @variable(deamodel, eff, start = 1.0)
        @variable(deamodel, lambda[1:nref] >= 0)

        @NLobjective(deamodel, Min, eff)

        @NLconstraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) <= eff^alpha * x0[j])
        @NLconstraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) >= y0[j] / (eff^(1-alpha)) )

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

    return GeneralizedDFDEAModel(n, m, s, alpha, rts, efficiency, lambdaeff)

end

function deagdf(X::Vector, Y::Matrix, alpha::Float64; rts::Symbol = :CRS, Xref::Vector = X, Yref::Matrix = Y)::GeneralizedDFDEAModel
    X = X[:,:]
    Xref = X[:,:]
    return deagdf(X, Y, alpha, rts = rts, Xref = Xref, Yref = Yref)
end

function deagdf(X::Matrix, Y::Vector, alpha::Float64; rts::Symbol = :CRS, Xref::Matrix = X, Yref::Vector = Y)::GeneralizedDFDEAModel
    Y = Y[:,:]
    Yref = Y[:,:]
    return deagdf(X, Y, alpha, rts = rts, Xref = Xref, Yref = Yref)
end

function deagdf(X::Vector, Y::Vector, alpha::Float64; rts::Symbol = :CRS, Xref::Vector = X, Yref::Vector = Y)::GeneralizedDFDEAModel
    X = X[:,:]
    Xref = X[:,:]
    Y = Y[:,:]
    Yref = Y[:,:]
    return deagdf(X, Y, alpha, rts = rts, Xref = Xref, Yref = Yref)
end

function Base.show(io::IO, x::GeneralizedDFDEAModel)
    compact = get(io, :compact, false)

    eff = efficiency(x)
    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)

    if !compact
        print(io, "Generalized DF DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "alpha = ", x.alpha)
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        show(io, CoefTable(hcat(eff), ["efficiency"], ["$i" for i in 1:n]))
    else

    end
end
