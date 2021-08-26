"""
    MalmquistDEAModel
An data structure representing a malmquist productivity DEA model.
"""
struct MalmquistDEAModel <: AbstractProductivityDEAModel
    n::Int64
    m::Int64
    s::Int64
    periods::Int64
    orient::Symbol
    rts::Symbol
    refperiod::Symbol
    dmunames::Union{Vector{String},Nothing}
    Prod::Matrix
    EC::Matrix
    TC::Matrix
end

"""
    malmquist(X, Y)
Compute the Malmquist productivity index using data envelopment analysis for inputs X and outputs Y.

# Optional Arguments
- `orient=:Input`: chooses between input oriented radial model `:Input` or output oriented radial model `:Output`.
- `refperiod=:Geomean`: chooses reference period for technological change: `:Base`, `:Comparison` or `:Geomean`.
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `names`: a vector of strings with the names of the decision making units.

# Examples
```jldoctest
julia> X = Array{Float64,3}(undef, 5, 1, 2);

julia> X[:, :, 1] = [2; 3; 5; 4; 4];

julia> X[:, :, 2] = [1; 2; 4; 3; 4];

julia> Y = Array{Float64,3}(undef, 5, 1, 2);

julia> Y[:, :, 1] = [1; 4; 6; 3; 5];

julia> Y[:, :, 2] = [1; 4; 6; 3; 3];

julia> malmquist(X, Y)
Mamlmquist DEA Model 
DMUs = 5; Inputs = 1; Outputs = 1; Time periods = 2
Orientation = Input; Returns to Scale = CRS
Referene period = Geomean
─────────────────────────
         M        EC   TC
─────────────────────────
1  2.0      1.33333   1.5
2  1.5      1.0       1.5
3  1.25     0.833333  1.5
4  1.33333  0.888889  1.5
5  0.6      0.4       1.5
─────────────────────────
M  = Malmquist Productivity Index 
EC = Efficiency Change 
TC = Technological Change
```
"""
function malmquist(X::Array{Float64,3}, Y::Array{Float64,3};
    orient::Symbol = :Input, rts::Symbol = :CRS, refperiod::Symbol = :Geomean,
    names::Union{Vector{String},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing)::MalmquistDEAModel

    # Check parameters
    nx, m, Tx = size(X)
    ny, s, Ty = size(Y)

    if nx != ny
        throw(DimensionMismatch("number of rows in X and Y ($nx, $ny) are not equal"));
    end
    if Tx != Ty
        throw(DimensionMismatch("number of time periods in X and Y ($Tx, $Ty) are not equal"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(GLPK.Optimizer)
    end

    # Compute efficiency for each DMU
    n = nx
    T = Tx

    Prod  = zeros(n, T - 1)
    EC    = zeros(n, T - 1)
    TC    = zeros(n, T - 1)

    # For each time period
    @showprogress 1 "Computing..." for t = 1:T - 1
        sleep(0.1)

        # Compute efficiency in t
        efft = efficiency(dea(X[:, :, t], Y[:, :, t], orient = orient, rts = rts, slack = false, optimizer = optimizer))

        # Compute efficiency in t + 1
        efft1 = efficiency(dea(X[:, :, t + 1], Y[:, :, t + 1], orient = orient, rts = rts, slack = false, optimizer = optimizer))

        # Compute efficiency in t +1, with reference t
        efft1_reft = efficiency(dea(X[:, :, t + 1], Y[:, :, t + 1], orient = orient, rts = rts, slack = false,
                                    Xref = X[:, :, t], Yref = Y[:, :, t], 
                                    optimizer = optimizer))

        # Compute efficiency in t, with reference t + 1
        efft_reft1 = efficiency(dea(X[:, :, t], Y[:, :, t], orient = orient, rts = rts, slack = false,
                                    Xref = X[:, :, t + 1], Yref = Y[:, :, t + 1], 
                                    optimizer = optimizer))

        # Inver if output oriented
        if orient == :Output
            efft       .= 1 ./ efft
            efft1      .= 1 ./ efft1
            efft1_reft .= 1 ./ efft1_reft
            efft_reft1 .= 1 ./ efft_reft1
        end

        # Technical Efficiency Change
        EC[:, t] .= efft1 ./ efft

        # Technological Change
        if refperiod == :Base
            TC[:, t] .= efft1_reft ./ efft1
        elseif refperiod == :Comparison
            TC[:, t] .= efft ./ efft_reft1
        elseif refperiod == :Geomean
            TC[:, t] .= sqrt.( (efft ./ efft_reft1) .* (efft1_reft ./ efft1) )
        else
            throw(ArgumentError("`refperiod` must be :Base, :Comparison or :Geomean"));
        end

        # Mamlmquist Index
        Prod[:, t] .= EC[:, t] .* TC[:, t]

    end

    return MalmquistDEAModel(n, m, s, T, orient, rts, refperiod, names, Prod, EC, TC)

end

function Base.show(io::IO, x::MalmquistDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    dmunames = names(x)

    periods = nperiods(x)
    Prod = prodchange(x, :Prod)
    EC   = prodchange(x, :EC)
    TC   = prodchange(x, :TC)

    if !compact
        print(io, "Mamlmquist DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "; Time periods = ", periods)
        print(io, "\n")
        print(io, "Orientation = ", string(x.orient))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        print(io, "Referene period = ", string(x.refperiod))
        print(io, "\n")
        if periods == 2
            show(io, CoefTable(hcat(Prod, EC, TC), ["M", "EC", "TC"], dmunames))
        end
        if periods > 2
            show(io, CoefTable(hcat(Prod, EC, TC),
                               vcat(["M$i" for i in 1:2], ["EC$i" for i in 1:2], ["TC$i" for i in 1:2 ]),
                               dmunames))
        end
        print(io, "\n")
        print(io, "M  = Malmquist Productivity Index \n")
        print(io, "EC = Efficiency Change \n")
        print(io, "TC = Technological Change")
    end

end
