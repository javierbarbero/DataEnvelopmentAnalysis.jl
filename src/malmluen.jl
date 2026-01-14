"""
    MalmquistLuenbergerDEAModel
An data structure representing a Malmquist-Luenberger productivity DEA model.
"""
struct MalmquistLuenbergerDEAModel <: AbstractProductivityDEAModel
    n::Int64
    m::Int64
    s::Int64
    b::Int64
    periods::Int64
    Gx::Symbol
    Gy::Symbol
    Gb::Symbol
    rts::Symbol
    refperiod::Symbol
    dmunames::Union{Vector{AbstractString},Nothing}
    Prod::Matrix
    EC::Matrix
    TC::Matrix
end

"""
    malmluen(X, Y, B; Gx, Gy, Gb)
Compute the Malmquist-Luenberger productivity index using data envelopment analysis for inputs `X`, 
good outputs `Y`, and bad outputs `B`, using directions `Gx`, `Gy`, and `Gb`.

# Direction specification:

The directions `Gx`, `Gy`, and `Gb` can be one of the following symbols.
- `:Zeros`: use zeros.
- `:Ones`: use ones.
- `:Observed`: use observed values.
- `:Mean`: use column means.

# Optional Arguments
- `refperiod=:Geomean`: chooses reference period for technological change: `:Base`, `:Comparison` or `:Geomean`.
- `names`: a vector of strings with the names of the decision making units.
"""
function malmluen(X::Array{Float64,3}, Y::Array{Float64,3}, B::Array{Float64,3};
    Gx::Union{Symbol, Array} = :Zeros, Gy::Union{Symbol, Array} = :Observed, Gb::Union{Symbol, Array} = :Observed,
    refperiod::Symbol = :Geomean,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::MalmquistLuenbergerDEAModel

    # Check parameters
    nx, m, Tx = size(X)
    ny, s, Ty = size(Y)
    nb, b, Tb = size(B)

    nx == ny || throw(DimensionMismatch("number of rows in X and Y ($nx, $ny) are not equal"));
    ny == nb || throw(DimensionMismatch("number of rows in Y and B ($ny, $nb) are not equal"));
    Tx == Ty || throw(DimensionMismatch("number of time periods in X and Y ($Tx, $Ty) are not equal"));
    Ty == Tb || throw(DimensionMismatch("number of time periods in Y and B ($Ty, $Tb) are not equal"));
    
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

    if typeof(Gb) == Symbol
        Gbsym = Gb

        if Gb == :Zeros
            Gb = zeros(size(B))
        elseif Gb == :Ones
            Gb = ones(size(B))
        elseif Gb == :Observed
            Gb = B
        elseif Gb == :Mean
            Gb = repeat(mean(B, dims = 1), size(B, 1))
        else
            throw(ArgumentError("Invalid `Gb`"));
        end
    else
        Gbsym = :Custom
    end
    
    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:LP)
    end

    # Compute efficiency for each DMU
    n = nx
    T = Tx

    Prod  = zeros(n, T - 1)
    EC    = zeros(n, T - 1)
    TC    = zeros(n, T - 1)

    # For each time period
    for t = 1:T - 1

        # Compute efficiency in t
        efft = efficiency(deaenv(X[:, :, t], Y[:, :, t], B[:, :, t], 
                                 Gx = Gx[:, :, t], Gy = Gy[:, :, t], Gb = Gb[:, :, t],
                                 slack = false, optimizer = optimizer))

        # Compute efficiency in t + 1
        efft1 = efficiency(deaenv(X[:, :, t + 1], Y[:, :, t + 1], B[:, :, t + 1], 
                                  Gx = Gx[:, :, t + 1], Gy = Gy[:, :, t + 1], Gb = Gb[:, :, t + 1],
                                  slack = false, optimizer = optimizer))

        # Compute efficiency in t +1, with reference t
        efft1_reft = efficiency(deaenv(X[:, :, t + 1], Y[:, :, t + 1], B[:, :, t + 1], 
                                    Gx = Gx[:, :, t + 1], Gy = Gy[:, :, t + 1], Gb = Gb[:, :, t + 1], 
                                    slack = false,
                                    Xref = X[:, :, t], Yref = Y[:, :, t], Bref = B[:, :, t], 
                                    optimizer = optimizer))

        # Compute efficiency in t, with reference t + 1
        efft_reft1 = efficiency(deaenv(X[:, :, t], Y[:, :, t], B[:, :, t],  
                                    Gx = Gx[:, :, t], Gy = Gy[:, :, t], Gb = Gb[:, :, t], 
                                    slack = false,
                                    Xref = X[:, :, t + 1], Yref = Y[:, :, t + 1], Bref = B[:, :, t + 1], 
                                    optimizer = optimizer))

        # Technical Efficiency Change
        EC[:, t] .= (1 .+ efft) ./ (1 .+ efft1)

        # Technological Change
        if refperiod == :Base
            TC[:, t] .= (1 .+ efft1) ./ (1 .+ efft1_reft)
        elseif refperiod == :Comparison
            TC[:, t] .= (1 .+ efft_reft1) ./ (1 .+ efft)
        elseif refperiod == :Geomean
            TC[:, t] .= sqrt.( ((1 .+ efft1) ./ (1 .+ efft1_reft)) .* ((1 .+ efft_reft1) ./ (1 .+ efft)) )
        else
            throw(ArgumentError("`refperiod` must be :Base, :Comparison or :Geomean"));
        end

        # Mamlmquist Index
        Prod[:, t] .= EC[:, t] .* TC[:, t]

    end

    return MalmquistLuenbergerDEAModel(n, m, s, b, T, Gxsym, Gysym, Gbsym, :CRS, refperiod, names, Prod, EC, TC)

end

isenvironmental(::MalmquistLuenbergerDEAModel) = true;

function Base.show(io::IO, x::MalmquistLuenbergerDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    b = nbadoutputs(x)
    dmunames = names(x)

    periods = nperiods(x)
    Prod = prodchange(x, :Prod)
    EC   = prodchange(x, :EC)
    TC   = prodchange(x, :TC)

    if !compact
        print(io, "Mamlmquist-Luenberger environmental DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "; Bad Outputs = ", b)
        print(io, "; Time periods = ", periods)
        print(io, "\n")
        print(io, "Returns to Scale = ", string(x.rts))
        print(io, "\n")
        print(io, "Gx = ", string(x.Gx), "; Gy = ", string(x.Gy), " ; Gb = ", string(x.Gb))
        print(io, "\n")
        print(io, "Referene period = ", string(x.refperiod))
        print(io, "\n")
        if periods == 2
            show(io, MIME"text/plain"(), CoefTable(hcat(Prod, EC, TC), ["ML", "EC", "TC"], dmunames))
        end
        if periods > 2
            show(io, MIME"text/plain"(), CoefTable(hcat(Prod, EC, TC),
                                         vcat(["ML$i" for i in 1:2], ["EC$i" for i in 1:2], ["TC$i" for i in 1:2 ]),
                                         dmunames))
        end
        print(io, "\n")
        print(io, "L  = Malmquist-Luenberger Productivity Index \n")
        print(io, "EC = Efficiency Change \n")
        print(io, "TC = Technological Change")
    end

end
