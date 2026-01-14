# This file contains functions for the Profit Efficiency DEA model
"""
    ProfitDEAModel
An data structure representing a profit DEA model.
"""
struct ProfitDEAModel <: AbstractProfitDEAModel
    n::Int64
    m::Int64
    s::Int64
    Gx::Symbol
    Gy::Symbol
    monetary::Bool
    dmunames::Union{Vector{AbstractString},Nothing}
    eff::Vector
    lambda::SparseMatrixCSC{Float64, Int64}
    techeff::Vector
    alloceff::Vector
    normalization::Vector
    Xtarget::Matrix
    Ytarget::Matrix
end


"""
    deaprofit(X, Y, W, P; Gx, Gy)
Compute profit efficiency using data envelopment analysis model for
inputs `X`, outputs `Y`, price of inputs `W`, and price of outputs `P`.

# Direction specification:

The directions `Gx` and `Gy` can be one of the following symbols.
- `:Zeros`: use zeros.
- `:Ones`: use ones.
- `:Observed`: use observed values.
- `:Mean`: use column means.
- `:Monetary`: use direction so that profit inefficiency is expressed in monetary values.
- `:Euclidean`: use prices normalized by the Euclidean norm of all input and output prices so profit inefficiency is the Euclidean distance between profit hyperplanes.

Alternatively, a vector or matrix with the desired directions can be supplied.

# Optional Arguments
- `monetary=false`: decomposition in normalized terms. Monetary terms if `true`.
- `names`: a vector of strings with the names of the decision making units.
"""
function deaprofit(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector},
    W::Union{Matrix,Vector}, P::Union{Matrix,Vector};
    Gx::Union{Symbol,Matrix,Vector}, Gy::Union{Symbol,Matrix,Vector},
    rts::Symbol = :VRS, monetary::Bool = false,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::ProfitDEAModel

    # Check parameters
    nx, m = size(X, 1), size(X, 2)
    ny, s = size(Y, 1), size(Y, 2)

    nw, mw = size(W, 1), size(W, 2)
    np, sp = size(P, 1), size(P, 2)

    if nx != ny
        throw(DimensionMismatch("number of rows in X and Y ($nx, $ny) are not equal"));
    end
    if nw != nx
        throw(DimensionMismatch("number of rows in W and X ($nw, $nx) are not equal"));
    end
    if np != ny
        throw(DimensionMismatch("number of rows in P and Y ($np, $ny) are not equal"));
    end
    if mw != m
        throw(DimensionMismatch("number of columns in W and X ($mw, $m) are not equal"));
    end
    if sp != s
        throw(DimensionMismatch("number of columns in P and Y ($sp, $s) are not equal"));
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
        elseif Gx == :Monetary
            GxGydollar = 1 ./ (sum(P, dims = 2) + sum(W, dims = 2));
            Gx = repeat(GxGydollar, 1, m);
        elseif Gx == :Euclidean
            Gx = W ./ sqrt.(sum(P.^2, dims = 2) + sum(W.^2, dims = 2));
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
        elseif Gy == :Monetary
            GxGydollar = 1 ./ (sum(P, dims = 2) + sum(W, dims = 2));
            Gy = repeat(GxGydollar, 1, s);
        elseif Gy == :Euclidean
            Gy = P ./ sqrt.(sum(P.^2, dims = 2) + sum(W.^2, dims = 2));
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

    # Get maximum profit targets and lambdas
    n = nx

    Xtarget, Ytarget, plambdaeff = deamaxprofit(X, Y, W, P, rts = rts, optimizer = optimizer)

    # Profit, technical and allocative efficiency
    maxprofit = sum(P .* Ytarget, dims = 2) .- sum(W .* Xtarget, dims = 2)

    pefficiency  = vec(maxprofit .- ( sum(P .* Y, dims = 2) .- sum(W .* X, dims = 2)))
    normalization = vec(sum(P .* Gy, dims = 2) .+ sum(W .* Gx, dims = 2))
    techefficiency = efficiency(deaddf(X, Y, Gx = Gx, Gy = Gy, rts = rts, slack = false, optimizer = optimizer))

    if monetary
        techefficiency = techefficiency .* normalization
    else
        pefficiency = pefficiency ./ normalization
    end

    allocefficiency = pefficiency - techefficiency

    return ProfitDEAModel(n, m, s, Gxsym, Gysym, monetary, names, pefficiency, plambdaeff, techefficiency, allocefficiency, normalization, Xtarget, Ytarget)

end

function Base.show(io::IO, x::ProfitDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    dmunames = names(x)

    eff = efficiency(x)
    techeff = efficiency(x, :Technical)
    alloceff = efficiency(x, :Allocative)

    if !compact
        print(io, "Profit DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Returns to Scale = VRS")
        print(io, "\n")
        print(io, "Gx = ", string(x.Gx), "; Gy = ", string(x.Gy))
        print(io, "\n")
        show(io, MIME"text/plain"(), CoefTable(hcat(eff, techeff, alloceff), ["Profit", "Technical", "Allocative"], dmunames))
    end

end
