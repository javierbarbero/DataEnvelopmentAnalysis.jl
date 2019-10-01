# This file contains functions for the Additive DEA model
"""
    AdditivelDEAModel
An data structure representing an additive DEA model.
"""
struct AdditiveDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64
    rts::Symbol
    eff::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    weights::Symbol
end

"""
    deaadd(X, Y)
Compute data envelopment analysis weighted additive model for inputs `X` and outputs `Y`.

# Optional Arguments
- `rts=:CRS`: chooses variable returns to scale. For constant returns to scale choose `:CRS`.
- `wX=ones(size(X))`: matrix of weights of inputs.
- `wY=ones(size(Y))`: matrix of weights of outputs.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaadd(X, Y)
Weighted Additive DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Weights = Ones; Returns to Scale = VRS
────────────────────────────────────────────────────
      efficiency       slackX1  slackX2      slackY1
────────────────────────────────────────────────────
1    0.0           0.0              0.0  0.0
2    7.33333       4.33333          0.0  3.0
3    0.0           0.0              0.0  0.0
4   -8.03397e-16  -8.03397e-16      0.0  0.0
5   18.0          13.0              1.0  4.0
6    6.48305e-16   2.70127e-16      0.0  3.78178e-16
7    0.0           0.0              0.0  0.0
8    0.0           0.0              0.0  0.0
9    0.0           0.0              0.0  0.0
10  35.0          25.0             10.0  0.0
11   4.0           0.0              4.0  4.78849e-16
────────────────────────────────────────────────────
```
"""
function deaadd(X::Matrix, Y::Matrix; rts::Symbol = :VRS, wX::Matrix = ones(size(X)), wY::Matrix = ones(size(Y)), Xref::Matrix = X, Yref::Matrix = Y)::AdditiveDEAModel
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
    if size(wX) != size(X)
        error("size of weights matrix for inputs should be equal to size of inputs")
    end
    if size(wY) != size(Y)
        error("size of weights matrix for outputs should be qual to size of outputs")
    end

    # Compute efficiency for each DMU
    n = nx
    nref = nrefx

    effi = zeros(n)
    slackX = zeros(n, m)
    slackY = zeros(n, s)
    lambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Value of weights to evaluate
        wX0 = wX[i,:]
        wY0 = wY[i,:]

        # Create the optimization model
        deamodel = Model(with_optimizer(GLPK.Optimizer))
        @variable(deamodel, sX[1:m] >= 0)
        @variable(deamodel, sY[1:s] >= 0)
        @variable(deamodel, lambda[1:nref] >= 0)

        @objective(deamodel, Max, sum(wX0[j] * sX[j] for j in 1:m) + sum(wY0[j] * sY[j] for j in 1:s) )
        @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) == x0[j] - sX[j])
        @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) == y0[j] + sY[j])

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
        slackX[i,:] = JuMP.value.(sX)
        slackY[i,:] = JuMP.value.(sY)

    end

    if all(wX .== 1) & all(wY .== 1)
        weights = :Ones
    else
        weights = :Custom
    end

    return AdditiveDEAModel(n, m, s, rts, effi, slackX, slackY, lambdaeff, weights)

end

function deaadd(X::Vector, Y::Matrix; rts::Symbol = :VRS, wX::Vector = ones(size(X)), wY::Matrix = ones(size(Y)), Xref::Vector = X, Yref::Matrix = Y)::AdditiveDEAModel
    X = X[:,:]
    wX = wX[:,:]
    Xref = Xref[:,:]
    return deaadd(X, Y, rts = rts, wX = wX, wY = wY, Xref = Xref, Yref = Yref)
end


function deaadd(X::Matrix, Y::Vector; rts::Symbol = :VRS, wX::Matrix = ones(size(X)), wY::Vector = ones(size(Y)), Xref::Matrix = X, Yref::Vector = Y)::AdditiveDEAModel
    Y = Y[:,:]
    wY = wY[:,:]
    Yref = Yref[:,:]
    return deaadd(X, Y, rts = rts, wX = wX, wY = wY, Xref = Xref, Yref = Yref)
end

function deaadd(X::Vector, Y::Vector; rts::Symbol = :VRS, wX::Vector = ones(size(X)), wY::Vector = ones(size(Y)), Xref::Vector = X, Yref::Vector = Y)::AdditiveDEAModel
    X = X[:,:]
    wX = wX[:,:]
    Xref = Xref[:,:]
    Y = Y[:,:]
    wY = wY[:,:]
    Yref = Yref[:,:]
    return deaadd(X, Y, rts = rts, wX = wX, wY = wY, Xref = Xref, Yref = Yref)
end

"""
    deaadd(X, Y. model)
Compute related data envelopment analysis weighted additive models for inputs `X` and outputs `Y`.

Model specification:
- `:Ones`: standard additive DEA model.
- `:MIP`: Measure of Inefficiency Proportions. (Charnes et al., 1987; Cooper et al., 1999)
- `:Normalized`: Normalized weighted additive DEA model. (Lovell and Pastor, 1995)
- `:RAM`: Range Adjusted Measure. (Cooper et al., 1999)
- `:BAM`: Bounded Adjusted Measure. (Cooper et al, 2011)

# Optional Arguments
- `rts=:VRS`: chosse between constant returns to scale `:CRS` or variable returns to scale `:VRS`.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaadd(X, Y, :MIP)
Weighted Additive DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Weights = MIP; Returns to Scale = VRS
─────────────────────────────────────────────────────
      efficiency       slackX1  slackX2       slackY1
─────────────────────────────────────────────────────
1    0.0           0.0              0.0   0.0
2    0.507519      0.0              0.0   7.10526
3    0.0           0.0              0.0   0.0
4   -4.72586e-17  -8.03397e-16      0.0   0.0
5    2.20395       0.0              0.0  17.6316
6    1.31279e-16   8.10382e-16      0.0   8.64407e-16
7    0.0           0.0              0.0   0.0
8    0.0           0.0              0.0   0.0
9    0.0           0.0              0.0   0.0
10   1.04322      17.0             15.0   1.0
11   0.235294      0.0              4.0   0.0
─────────────────────────────────────────────────────
```
"""
function deaadd(X::Matrix, Y::Matrix, model::Symbol; rts::Symbol = :VRS, Xref::Matrix = X, Yref::Matrix = Y)::AdditiveDEAModel
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

    # Compute specific weights based on the model
    if model == :Ones
        # Standard Additive DEA model
        wX = ones(size(X))
        wY = ones(size(Y))
    elseif model == :MIP
        # Measure of Inefficiency Proportions
        wX = 1 ./ X
        wY = 1 ./ Y
    elseif model == :Normalized
        # Normalized weighted additive DEA model
        wX = zeros(size(X))
        wY = zeros(size(Y))

        m = size(X, 2)
        s = size(Y, 2)

        for i=1:m
            wX[:,i] .= 1 ./ std(Xref[:,i])
        end
        for i=1:s
            wY[:,i] .= 1 ./ std(Yref[:,i])
        end

    elseif model == :RAM
        # Range Adjusted Measure
        m = size(X, 2)
        s = size(Y, 2)

        wX = zeros(size(X))
        wY = zeros(size(Y))

        for i=1:m
            wX[:,i] .= 1 ./ ((m + s) * (maximum(Xref[:,i])  - minimum(Xref[:,i])))
        end
        for i=1:s
            wY[:,i] .= 1 ./ ((m + s) * (maximum(Yref[:,i])  - minimum(Yref[:,i])))
        end

    elseif model == :BAM
        # Bounded Adjusted Measure
        m = size(X, 2)
        s = size(Y, 2)

        minXref = zeros(m)
        maxYref = zeros(s)

        wX = zeros(size(X))
        wY = zeros(size(Y))

        for i=1:m
            minXref[i] = minimum(Xref[:,i])
            wX[:,i] = 1 ./ ((m  + s) .* (X[:,i] .- minXref[i] ))
        end
        for i=1:s
            maxYref[i] = maximum(Yref[:,i])
            wY[:,i] = 1 ./ ((m + s) .* (maxYref[i] .- Y[:,i]))
        end

        wX[isinf.(wX)] .= 0
        wY[isinf.(wY)] .= 0

    else
        error("Invalid model ", model)
    end

    # Compute efficiency for each DMU
    n = nx
    nref = nrefx

    effi = zeros(n)
    slackX = zeros(n, m)
    slackY = zeros(n, s)
    lambdaeff = spzeros(n, nref)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Value of weights to evaluate
        wX0 = wX[i,:]
        wY0 = wY[i,:]

        # Create the optimization model
        deamodel = Model(with_optimizer(GLPK.Optimizer))
        @variable(deamodel, sX[1:m] >= 0)
        @variable(deamodel, sY[1:s] >= 0)
        @variable(deamodel, lambda[1:nref] >= 0)

        @objective(deamodel, Max, sum(wX0[j] * sX[j] for j in 1:m) + sum(wY0[j] * sY[j] for j in 1:s) )
        @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) == x0[j] - sX[j])
        @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) == y0[j] + sY[j])

        # Add return to scale constraints
        if rts == :CRS
            # Add constraints for BAM CRS model
            if model == :BAM
                @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * lambda[t] for t in 1:nref) >= minXref[j])
                @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * lambda[t] for t in 1:nref) <= maxYref[j])
            end
        elseif rts == :VRS
            @constraint(deamodel, sum(lambda) == 1)
        else
            error("Invalid returns to scale $rts. Returns to scale should be :CRS or :VRS")
        end

        # Optimize and return results
        JuMP.optimize!(deamodel)

        effi[i]  = JuMP.objective_value(deamodel)
        lambdaeff[i,:] = JuMP.value.(lambda)
        slackX[i,:] = JuMP.value.(sX)
        slackY[i,:] = JuMP.value.(sY)

    end

    return AdditiveDEAModel(n, m, s, rts, effi, slackX, slackY, lambdaeff, model)

end

function deaadd(X::Vector, Y::Matrix, model::Symbol; rts::Symbol = :VRS, Xref::Vector = X, Yref::Matrix = Y)::AdditiveDEAModel
    X = X[:,:]
    Xref = Xref[:,:]
    return deaadd(X, Y, model, rts = rts, Xref = Xref, Yref = Yref)
end

function deaadd(X::Matrix, Y::Vector, model::Symbol; rts::Symbol = :VRS, Xref::Matrix = X, Yref::Vector = Y)::AdditiveDEAModel
    Y = Y[:,:]
    Yref = Yref[:,:]
    return deaadd(X, Y, model, rts = rts, Xref = Xref, Yref = Yref)
end

function deaadd(X::Vector, Y::Vector, model::Symbol; rts::Symbol = :VRS, Xref::Vector = X, Yref::Vector = Y)::AdditiveDEAModel
    X = X[:,:]
    Xref = Xref[:,:]
    Y = Y[:,:]
    Yref = Yref[:,:]
    return deaadd(X, Y, model, rts = rts, Xref = Xref, Yref = Yref)
end

function Base.show(io::IO, x::AdditiveDEAModel)
    compact = get(io, :compact, false)
    
    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    eff = efficiency(x)
    slackX = slacks(x, :X)
    slackY = slacks(x, :Y)

    if !compact
        print(io, "Weighted Additive DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Weights = ", string(x.weights))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        show(io, CoefTable(hcat(eff, slackX, slackY), ["efficiency"; ["slackX$i" for i in 1:m ]; ; ["slackY$i" for i in 1:s ]], ["$i" for i in 1:n]))
    end
    
end
