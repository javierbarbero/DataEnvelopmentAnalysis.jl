# This file contains functions for the Reverse DDF DEA model
"""
    ReverseDDFDEAModel
An data structure representing a Reverse DDF DEA model.
"""
struct ReverseDDFDEAModel <: AbstractTechnicalDEAModel
    n::Int64
    m::Int64
    s::Int64    
    measure::Symbol
    orient::Symbol
    rts::Symbol    
    dmunames::Union{Vector{String},Nothing}
    eff::Vector
    slackX::Matrix
    slackY::Matrix
    lambda::SparseMatrixCSC{Float64, Int64}
    Xtarget::Matrix
    Ytarget::Matrix
    Gxrddf::Union{Vector,Matrix}
    Gyrddf::Union{Vector,Matrix}
end

"""
    dearddf(X, Y, measure)
Compute data envelopment analysis reverse directional distance function (RDDF) model for inputs
`X`, outputs `Y`, and efficiency measure  `measure`.

# Measure specification:

- `:ERG`: Enhanced Russell Graph (or Slack Based Measure (SBM)).
- `:MDDF`: Modified Directional Distance Function.

# Direction specification:

For the Modified Directional Distance Function, the directions `Gx` and `Gy` can be one of the following symbols.
- `:Ones`: use ones.
- `:Observed`: use observed values.
- `:Mean`: use column means.

# Optional Arguments
- `orient=:Graph`: choose between graph oriented `:Graph`, input oriented `:Input`, or output oriented model `:Output`.
- `rts=:CRS`: choose between constant returns to scale `:CRS` or variable returns to scale `:VRS`.
- `atol=1e-6`: tolerance for DMU to be considered efficient.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `names`: a vector of strings with the names of the decision making units.
"""
function dearddf(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector}, measure::Symbol;
    Gx::Union{Symbol,Matrix,Vector,Nothing} = nothing, Gy::Union{Symbol,Matrix,Vector,Nothing} = nothing,
    orient::Symbol = :Graph, rts::Symbol = :CRS, slack::Bool = false, atol::Float64 = 1e-6, 
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix,Vector,Nothing} = nothing,
    names::Union{Vector{String},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::ReverseDDFDEAModel

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

    n = nx

    # Get directions based on efficiency measure
    Gxrddf, Gyrddf = dearddfdirections(X, Y, measure, Gx = Gx, Gy = Gy, orient = orient, rts = rts, atol = atol, Xref = Xref, Yref = Yref, optimizer = optimizer)

    # Directional model using calculated directions
    rddfmodel = deaddf(X, Y, Gx = Gxrddf, Gy = Gyrddf, rts = rts, Xref = Xref, Yref = Yref, slack = slack, optimizer = optimizer)
    
    effi = efficiency(rddfmodel)
    lambdaeff = peersmatrix(rddfmodel)
    slackX = slacks(rddfmodel, :X)
    slackY = slacks(rddfmodel, :Y)
    Xtarget = targets(rddfmodel, :X)
    Ytarget = targets(rddfmodel, :Y)

    return ReverseDDFDEAModel(n, m, s, measure, orient, rts, names, effi, slackX, slackY, lambdaeff, Xtarget, Ytarget, Gxrddf, Gyrddf)

end

"""
    dearddfdirections(X, Y, measure)
Compute corresponding directions for data envelopment analysis Reverse DDF models 
for inputs `X`, outputs `Y`, and efficiency `measure`.

# Measure specification:
- `:ERG`: Enhanced Russell Graph (or Slack Based Measure (SBM)).
- `:MDDF`: Modified Directional Distance Function.

# Direction specification:

For the Modified Directional Distance Function, the directions `Gx` and `Gy` can be one of the following symbols.
- `:Ones`: use ones.
- `:Observed`: use observed values.
- `:Mean`: use column means.

Alternatively, a vector or matrix with the desired directions can be supplied.
# Optional Arguments
- `orient=:Graph`: choose between graph oriented `:Graph`, input oriented `:Input`, or output oriented model `:Output`.
- `rts=:CRS`: choose between constant returns to scale `:CRS` or variable returns to scale `:VRS`.
- `atol=1e-6`: tolerance for DMU to be considered efficient.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
"""
function dearddfdirections(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector}, measure::Symbol;
    Gx::Union{Symbol,Matrix,Vector,Nothing} = nothing, Gy::Union{Symbol,Matrix,Vector,Nothing} = nothing,
    orient::Symbol = :Graph, rts::Symbol = :CRS, atol::Float64 = 1e-6,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix,Vector,Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)

    # Calculate direction based on efficiency measure
    Gxrddf = zeros(size(X))
    Gyrddf = zeros(size(Y))

    n = size(X, 1)

    if measure == :ERG
        # Show error if user specify directions
        if (Gx !== nothing) | (Gy !== nothing)
            throw(ArgumentError("`Gx` and `Gy` should not be specified for ERG measure."));
        end

        # Enhanced Russels Graph  (ERG) (or Slack Based Measure (SBM))
        if orient == :Graph
            ergmodel = deaerg(X, Y, rts = rts, Xref = Xref, Yref = Yref, optimizer = optimizer)
        elseif orient == :Input
            ergmodel = dearussell(X, Y, orient = :Input, rts = rts, Xref = Xref, Yref = Yref, slack = false, optimizer = optimizer)
        elseif orient == :Output
            ergmodel = dearussell(X, Y, orient = :Output, rts = rts, Xref = Xref, Yref = Yref, slack = false, optimizer = optimizer)
        else
            throw(ArgumentError("`orient` must be :Input or :Output"));
        end

        eff = efficiency(ergmodel)
        if orient == :Graph || orient == :Input
            ineff = 1 .- eff
        elseif orient == :Output
            #ineff = eff .- 1
            ineff = 1 .- (1 ./ eff)
        end

        Xtarget = targets(ergmodel, :X)
        Ytarget = targets(ergmodel, :Y)

        for i = 1:n
            if abs(ineff[i]) <= atol
                Gxrddf[i,:] .= 1
                Gyrddf[i,:] .= 1
            else
                Gxrddf[i,:] = (X[i,:] - Xtarget[i,:]) ./ ineff[i]
                Gyrddf[i,:] = (Ytarget[i,:] - Y[i,:]) ./ ineff[i]
            end
        end
    elseif measure == :MDDF
        # Modified DDF
        if orient != :Graph
            throw(ArgumentError("Measure MDDF only available for Graph orientation"))
        end

        # Build or get user directions
        if typeof(Gx) == Symbol
            if Gx == :Ones
                Gx = ones(size(X))
            elseif Gx == :Observed
                Gx = X
            elseif Gx == :Mean
                Gx = repeat(mean(X, dims = 1), size(X, 1))
            else
                throw(ArgumentError("Invalid `Gx`"));
            end
        end

        if typeof(Gy) == Symbol
            if Gy == :Ones
                Gy = ones(size(Y))
            elseif Gy == :Observed
                Gy = Y
            elseif Gy == :Mean
                Gy = repeat(mean(Y, dims = 1), size(Y, 1))
            else
                throw(ArgumentError("Invalid `Gy`"));
            end
        end

        mddfmodel = deamddf(X, Y, Gx = Gx, Gy = Gy, rts = rts, Xref = Xref, Yref = Yref, slack = false, optimizer = optimizer)

        betax = mddfmodel.betax
        betay = mddfmodel.betay

        for i = 1:n
            if betax[i] != 0
                Gxrddf[i,:] = (betax[i] / (betax[i] + betay[i])) .* Gx[i,:]
            else
                Gxrddf[i,:] .= 0
            end

            if betay[i] != 0
                Gyrddf[i,:] = (betay[i] / (betax[i] + betay[i])) .* Gy[i,:]
            else
                Gyrddf[i,:] .= 0
            end
        end

    else
        throw(ArgumentError("Invalid efficiency `measure`"));
    end

    # Set directions to 0 if Input or Output oriented
    if orient == :Input
        Gyrddf = zeros(size(Y))
    elseif orient == :Output
        Gxrddf = zeros(size(X))
    end

    return Gxrddf, Gyrddf

end

function Base.show(io::IO, x::ReverseDDFDEAModel)
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
        print(io, "Reverse DDF DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = ", string(x.orient))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        print(io, "Associated efficiency measure = ", string(x.measure))
        print(io, "\n")

        if hasslacks == true
            show(io, CoefTable(hcat(eff, slackX, slackY), ["efficiency"; ["slackX$i" for i in 1:m ]; ["slackY$i" for i in 1:s ]], dmunames))
        else
            show(io, CoefTable(hcat(eff), ["efficiency"], dmunames))
        end
    end

end
