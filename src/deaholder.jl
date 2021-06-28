# This file contains functions for the Hölder DEA model
"""
    AbstractHolderDEAModel
An abstract type representing a Holder DEA model.
"""
abstract type AbstractHolderDEAModel <: AbstractTechnicalDEAModel end

"""
    HolderL1DEAModel
An data structure representing a Höler L1 DEA model.
"""
struct HolderL1DEAModel <: AbstractHolderDEAModel
    n::Int64
    m::Int64
    s::Int64
    orient::Symbol
    rts::Symbol
    isweighted::Bool
    dmunames::Union{Vector{String},Nothing}
    eff::Vector
    effmin::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    Xtarget::Matrix
    Ytarget::Matrix
end

"""
    HolderL2DEAModel
An data structure representing a Höler L2 DEA model.
"""
struct HolderL2DEAModel <: AbstractHolderDEAModel
    n::Int64
    m::Int64
    s::Int64
    orient::Symbol
    rts::Symbol
    isweighted::Bool
    dmunames::Union{Vector{String},Nothing}
    eff::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    Xtarget::Matrix
    Ytarget::Matrix
end

"""
    HolderLInfDEAModel
An data structure representing a Höler LInf DEA model.
"""
struct HolderLInfDEAModel <: AbstractHolderDEAModel
    n::Int64
    m::Int64
    s::Int64
    orient::Symbol
    rts::Symbol
    isweighted::Bool
    dmunames::Union{Vector{String},Nothing}
    eff::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    Xtarget::Matrix
    Ytarget::Matrix
end

"""
    deaholder(X, Y; l)
Compute the Hölder distance function model using data envelopment analysis for inputs `X` and outputs `Y`, 
using Hölder norm `l`.

# Hölder norm `l` specification
- `1`.
- `2`.
- `Inf`.

# Optional Arguments
- `weigt=false`:  set to `true` for weighted (weakly) Hölder distance function.
- `orient=:Input`: chooses the radially oriented input mode. For the radially oriented output model choose `:Output`.
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `slack=true`: computes input and output slacks.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `names`: a vector of strings with the names of the decision making units.

# Examples
```jldoctest
julia> X = [2; 4; 8; 12; 6; 14; 14; 9.412];

julia> Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

julia> deaholder(X, Y, l = 1, rts = :VRS)
Hölder L1 DEA Model 
DMUs = 8; Inputs = 1; Outputs = 1
Orientation = Graph; Returns to Scale = VRS
────────────────────────────────────────────────
   efficiency  minimum      slackX1      slackY1
────────────────────────────────────────────────
1         0.0       X1  0.0          0.0
2         0.0       X1  0.0          0.0
3         0.0       X1  0.0          0.0
4         0.0       X1  0.0          0.0
5         3.0       X1  0.0          1.01506e-15
6         2.0       Y1  2.0          0.0
7         0.0       Y1  2.0          0.0
8         6.0       Y1  1.77636e-15  0.0
────────────────────────────────────────────────
```
"""
function deaholder(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    l::Union{Int64,Float64}, weight::Bool = false,
    orient::Symbol = :Graph, rts::Symbol = :CRS, slack = true,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix, Vector,Nothing} = nothing,
    names::Union{Vector{String},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::AbstractHolderDEAModel

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

    if (orient != :Input) && (orient != :Output) && (orient != :Graph)
        throw(ArgumentError("`orient` must be :Input or :Output"));
    end
    
    # Check l
    if l == 1
        return deaholderl1(X, Y, weight = weight, orient = orient, rts = rts, slack = slack, Xref = Xref, Yref = Yref, names = names, optimizer = optimizer)
    elseif l == 2
        return deaholderl2(X, Y, weight = weight, orient = orient, rts = rts, slack = slack, Xref = Xref, Yref = Yref, names = names, optimizer = optimizer)
    elseif l == Inf
        return deaholderlinf(X, Y, weight = weight, orient = orient, rts = rts, slack = slack, Xref = Xref, Yref = Yref, names = names, optimizer = optimizer)
    else
        throw(ArgumentError("`l` must be 1, 2 or Inf"));
    end
    
end

function deaholderl1(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    weight::Bool = false,
    orient::Symbol = :Graph, rts::Symbol = :CRS, slack = true,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix, Vector,Nothing} = nothing,
    names::Union{Vector{String},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::HolderL1DEAModel
    
    # Get parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    if Xref === nothing Xref = X end
    if Yref === nothing Yref = Y end

    nrefx, mref = size(Xref, 1), size(Xref, 2)
    nrefy, sref = size(Yref, 1), size(Yref, 2)

    # Compute DDF efficiency with 1 for each input and output
    n = nx
    nref = nrefx

    effXY = fill(Inf, (n,m+s))
    lambdaXY = Array{SparseMatrixCSC{Float64, Int64}}(undef, m + s)

    if slack == true
        slackXms = Array{Matrix}(undef, m + s)
        slackYms = Array{Matrix}(undef, m + s)
    end

    if orient == :Graph || orient == :Input
        for i = 1:m
            Gxi = zeros(n,m)
            if weight
                Gxi[:,i] .= X[:,i]
            else
                Gxi[:,i] .= 1
            end

            tempdea = deaddf(X, Y, Gx = Gxi, Gy = :Zeros, rts = rts, Xref = Xref, Yref = Yref, optimizer = optimizer)
            effXY[:,i] = efficiency(tempdea)
            lambdaXY[i] = tempdea.lambda

            if slack == true
                slackXms[i] = slacks(tempdea, :X)
                slackYms[i] = slacks(tempdea, :Y)
            end
        end
    end

    if orient == :Graph || orient == :Output
        for i = 1:s
            Gyi = zeros(n,s)
            if weight
                Gyi[:,i] .= Y[:,i]
            else
                Gyi[:,i] .= 1
            end

            tempdea = deaddf(X, Y, Gx = :Zeros, Gy = Gyi, rts = rts, Xref = Xref, Yref = Yref, optimizer = optimizer)
            effXY[:,i + m] = efficiency(tempdea)
            lambdaXY[i + m] = tempdea.lambda

            if slack == true
                slackXms[i + m] = slacks(tempdea, :X)
                slackYms[i + m] = slacks(tempdea, :Y)
            end
        end
    end

    # Obtain the minimum efficiency scores
    effi = zeros(n)
    effmin = fill(0, n)
    lambda = spzeros(n, nref)

    if slack == true
        slackX = zeros(n, m)
        slackY = zeros(n, s)
    else
        slackX = Array{Float64}(undef, 0, 0)
        slackY = Array{Float64}(undef, 0, 0)
    end

    for i = 1:n
        effi[i], effmin[i] = findmin(effXY[i,:])
        lambda[i, :] = lambdaXY[effmin[i]][i,:]

        if slack == true
            slackX[i,:] = slackXms[effmin[i]][i,:]
            slackY[i,:] = slackYms[effmin[i]][i,:]
        end
    end

    # Get X and Y targets
    Xtarget = convert.(Float64, X[:,:])
    Ytarget = convert.(Float64, Y[:,:])
    for i = 1:n
        if effmin[i] <= m
            if weight
                Xtarget[i,effmin[i]] = X[i,effmin[i]] - effi[i] * X[i,effmin[i]]
            else
                Xtarget[i,effmin[i]] = X[i,effmin[i]] - effi[i]
            end
        else
            if weight
                Ytarget[i,effmin[i] - m] = Y[i,effmin[i] - m] + effi[i] * Y[i,effmin[i] - m]
            else
                Ytarget[i,effmin[i] - m] = Y[i,effmin[i] - m] + effi[i]
            end
        end
    end

    if slack == true
        Xtarget = Xtarget - slackX
        Ytarget = Ytarget + slackY
    end
    
    return HolderL1DEAModel(n, m, s, orient, rts, weight, names,  effi, effmin, slackX, slackY, lambda, Xtarget, Ytarget)

end

function deaholderl2(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    weight::Bool = false,
    orient::Symbol = :Graph, rts::Symbol = :CRS, slack = true,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix, Vector,Nothing} = nothing,
    names::Union{Vector{String},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::HolderL2DEAModel

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    if Xref === nothing Xref = X end
    if Yref === nothing Yref = Y end

    nrefx, mref = size(Xref, 1), size(Xref, 2)
    nrefy, sref = size(Yref, 1), size(Yref, 2)

    if rts != :VRS
        throw(ArgumentError("`rts` must be :VRS if l = 2"));
    end

    # Default optimizer
    if optimizer === nothing         
        throw(ArgumentError("An optimizer that supports SOS constraints must be specified for Hölder l = 2"));
    end

    noSOS = false

    # Compute efficiency for each DMU
    n = nx
    nref = n

    effi = zeros(n)
    lambdaeff = spzeros(n, nref)
    Xtarget = convert.(Float64, X[:,:])
    Ytarget = convert.(Float64, Y[:,:])
    slackX = zeros(n, m)
    slackY = zeros(n, s)

    for i=1:n
        # Value of inputs and outputs to evaluate
        x0 = X[i,:]
        y0 = Y[i,:]

        # Create the optimization model
        deamodel = newdeamodel(optimizer)

        if orient == :Input || orient == :Graph
            @variable(deamodel, Xeff[1:m] >= 0)
            set_start_value.(Xeff, x0)
        end
        if orient == :Output || orient == :Graph
            @variable(deamodel, Yeff[1:s] >= 0)
            set_start_value.(Yeff, y0)
        end

        if weight
            if orient == :Graph
                @objective(deamodel, Min, sum(((x0[j] - Xeff[j]) / x0[j])^2 for j in 1:m) + sum(((Yeff[j] - y0[j]) / y0[j])^2 for j in 1:s)  )
            elseif orient == :Input
                @objective(deamodel, Min, sum(((x0[j] - Xeff[j]) / x0[j])^2 for j in 1:m) )
            elseif orient == :Output
                @objective(deamodel, Min, sum(((Yeff[j] - y0[j]) / y0[j])^2 for j in 1:s)  )
            end
        else
            if orient == :Graph
                @objective(deamodel, Min, sum((x0[j] - Xeff[j])^2 for j in 1:m) + sum((Yeff[j] - y0[j])^2 for j in 1:s)  )
            elseif orient == :Input
                @objective(deamodel, Min, sum((x0[j] - Xeff[j])^2 for j in 1:m) )
            elseif orient == :Output
                @objective(deamodel, Min, sum((Yeff[j] - y0[j])^2 for j in 1:s)  )
            end
        end
   
        @variable(deamodel, gammatau[1:2, 1:nref] >= 0)
        @variable(deamodel, vsX[1:2,1:m] >= 0)
        @variable(deamodel, musY[1:2,1:s] >= 0)
        @variable(deamodel, psi)

        set_start_value.(gammatau[1,i], 1)
        
        if orient == :Input || orient == :Graph
            @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * gammatau[1,t] for t in 1:nref) == Xeff[j] - vsX[2,j])
        elseif orient == :Output
            @constraint(deamodel, [j in 1:m], sum(Xref[t,j] * gammatau[1,t] for t in 1:nref) == x0[j] - vsX[2,j])
        end
        if orient == :Output || orient == :Graph
            @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * gammatau[1,t] for t in 1:nref) == Yeff[j] + musY[2,j])
        elseif orient == :Input
            @constraint(deamodel, [j in 1:s], sum(Yref[t,j] * gammatau[1,t] for t in 1:nref) == y0[j] + musY[2,j])
        end

        @constraint(deamodel, sum(gammatau[1,:]) == 1)

        if orient == :Graph
            @constraint(deamodel, sum(vsX[1,t] for t in 1:m) + sum(musY[1,t] for t in 1:s) == 1)
        elseif orient == :Input
            @constraint(deamodel, sum(vsX[1,t] for t in 1:m) == 1)
        elseif orient == :Output
            @constraint(deamodel, sum(musY[1,t] for t in 1:s) == 1)
        end

        @constraint(deamodel, [j in 1:nref], sum(vsX[1,t] * Xref[j,t] for t in 1:m) - sum(musY[1,t] * Yref[j,t] for t in 1:s) + psi - gammatau[2,j] == 0 )
     
        # SOS constraints
        try
            @constraint(deamodel, [j in 1:nref], gammatau[:,j] in SOS1())                           
        catch
            @constraint(deamodel, [j in 1:nref], gammatau[1,j] * gammatau[2,j] == 0)  
            noSOS = true
        end

        try
            @constraint(deamodel, [j in 1:m], vsX[:,j] in SOS1())
        catch
            @constraint(deamodel, [j in 1:m], vsX[1,j] * vsX[2,j] == 0)
            noSOS = true
        end

        try
            @constraint(deamodel, [j in 1:s], musY[:,j] in SOS1())    
        catch
            @constraint(deamodel, [j in 1:s], musY[1,j] * musY[2,j] == 0)   
            noSOS = true
        end
        
        # Optimize and return results
        JuMP.optimize!(deamodel)

        effi[i] = JuMP.objective_value(deamodel)
        if effi[i] < 0
            # If objective function is negative change to 0, to avoid error when calculatin the square root.
            @warn ("DMU $i negative value in the objective function: $(effi[i]). Replaced to 0")
            effi[i] = 0
        end
        effi[i] = sqrt(effi[i] )
        lambdaeff[i,:] = JuMP.value.(gammatau[1,:])

        if orient == :Input || orient == :Graph
            Xtarget[i,:]  = JuMP.value.(Xeff)
        end
        if orient == :Output || orient == :Graph
            Ytarget[i,:]  = JuMP.value.(Yeff)
        end

        slackX[i,:] = JuMP.value.(vsX[2,:])
        slackY[i,:] = JuMP.value.(musY[2,:])

        if (termination_status(deamodel) != MOI.OPTIMAL) && (termination_status(deamodel) != MOI.LOCALLY_SOLVED)
            @warn ("DMU $i termination status: $(termination_status(deamodel)). Primal status: $(primal_status(deamodel)). Dual status: $(dual_status(deamodel))")
        end

    end

    # Add slacks to X and Y targets
    Xtarget = Xtarget - slackX
    Ytarget = Ytarget + slackY

    # Display warning if model not solved with SOS constraints
    if noSOS
        @warn ("Model solved with $(optimizer.optimizer) is innacuate. Use a solver that supports SOS1 constraints.")
    end

    return HolderL2DEAModel(n, m, s, orient, rts, weight, names, effi, slackX, slackY, lambdaeff, Xtarget, Ytarget)

end


function deaholderlinf(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    weight::Bool = true,
    orient::Symbol = :Graph, rts::Symbol = :CRS, slack = true,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix, Vector,Nothing} = nothing,
    names::Union{Vector{String},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::HolderLInfDEAModel

    # Get parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    n = nx

    # Buil direction based on orientation
    if orient == :Input
        if weight
            Gx = X
        else
            Gx = :Ones
        end
        Gy = :Zeros
    elseif orient == :Output
        Gx = :Zeros
        if weight
            Gy = Y
        else
            Gy = :Ones
        end
    elseif orient == :Graph
        if weight
            Gx = X
            Gy = Y
        else
            Gx = :Ones
            Gy = :Ones
        end
    end

    # Get DDF technical efficiency
    tempdea = deaddf(X, Y, Gx = Gx, Gy = Gy, rts = rts, slack = slack, Xref = Xref, Yref = Yref, optimizer = optimizer)
    effi = efficiency(tempdea)
    lambda = peersmatrix(tempdea)

    slackX = slacks(tempdea, :X)
    slackY = slacks(tempdea, :Y)

    Xtarget = targets(tempdea, :X)
    Ytarget = targets(tempdea, :Y)

    return HolderLInfDEAModel(n, m, s, orient, rts, weight, names, effi, slackX, slackY, lambda, Xtarget, Ytarget)        

end

function efficiency(model::HolderL1DEAModel, type::Symbol)::Vector

    if (type == :min) return model.effmin end

    throw(ArgumentError("$(typeof(model)) with orienation has no efficiency $(type)"));

end

function Base.show(io::IO, x::HolderL1DEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    isweighted = x.isweighted
    dmunames = names(x)

    eff = efficiency(x)
    effmin = x.effmin
    slackX = slacks(x, :X)
    slackY = slacks(x, :Y)
    hasslacks = ! isempty(slackX)

    effmincol = Array{String}(undef,n)
    for i = 1:n
        if effmin[i] <= m
            effmincol[i] = "X" * string(Int(effmin[i]))
        else
            effmincol[i] = "Y" * string(Int(effmin[i] - m))
        end
    end

    if !compact
        print(io, "Hölder L1 DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = ", string(x.orient))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        if isweighted
            print(io, "Weighted (weakly) Hölder distance function \n")
        end

        if hasslacks == true
            show(io, CoefTable(hcat(eff, effmincol, slackX, slackY), ["efficiency"; "minimum"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
        else
            show(io, CoefTable(hcat(eff, effmincol), ["efficiency"; "minimum"], dmunames))
        end        
    end

end

function Base.show(io::IO, x::HolderL2DEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    isweighted = x.isweighted
    dmunames = names(x)

    eff = efficiency(x)
    slackX = slacks(x, :X)
    slackY = slacks(x, :Y)
    hasslacks = ! isempty(slackX)

    if !compact
        print(io, "Hölder L2 DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = ", string(x.orient))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        if isweighted
            print(io, "Weighted (weakly) Hölder distance function \n")
        end

        show(io, CoefTable(hcat(eff, slackX, slackY), ["efficiency"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
    end

end

function Base.show(io::IO, x::HolderLInfDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    isweighted = x.isweighted
    dmunames = names(x)

    eff = efficiency(x)
    slackX = slacks(x, :X)
    slackY = slacks(x, :Y)
    hasslacks = ! isempty(slackX)

    if !compact
        print(io, "Hölder LInf DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = ", string(x.orient))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        if isweighted
            print(io, "Weighted (weakly) Hölder distance function \n")
        end

        if hasslacks == true
            show(io, CoefTable(hcat(eff, slackX, slackY), ["efficiency"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
        else
            show(io, CoefTable(hcat(eff), ["efficiency"], dmunames))
        end
    end

end

