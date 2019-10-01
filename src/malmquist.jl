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
function malmquist(X::Array{Float64,3}, Y::Array{Float64,3}; orient::Symbol = :Input, rts::Symbol = :CRS, refperiod::Symbol = :Geomean)::MalmquistDEAModel
    # Check dimensions
    if ndims(X) != 3
        error("X should be a 3-dimensions array")
    end
    if ndims(Y) != 3
        error("Y should be a 3-dimensions array")
    end

    # Check parameters
    nx, m, Tx = size(X)
    ny, s, Ty = size(Y)

    if nx != ny
        error("number of observations is different in inputs and outputs")
    end
    if Tx != Ty
        error("number of time periods is different in intputs and outputs")
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
        efft = efficiency(dea(X[:, :, t], Y[:, :, t], orient = orient, rts = rts, slack = false))

        # Compute efficiency in t + 1
        efft1 = efficiency(dea(X[:, :, t + 1], Y[:, :, t + 1], orient = orient, rts = rts, slack = false))

        # Compute efficiency in t +1, with reference t
        efft1_reft = efficiency(dea(X[:, :, t + 1], Y[:, :, t + 1], orient = orient, rts = rts, slack = false,
                                    Xref = X[:, :, t], Yref = Y[:, :, t]))

        # Compute efficiency in t, with reference t + 1
        efft_reft1 = efficiency(dea(X[:, :, t], Y[:, :, t], orient = orient, rts = rts, slack = false,
                                    Xref = X[:, :, t + 1], Yref = Y[:, :, t + 1]))

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
        end

        # Mamlmquist Index
        Prod[:, t] .= EC[:, t] .* TC[:, t]

    end

    return MalmquistDEAModel(n, m, s, T, orient, rts, refperiod, Prod, EC, TC)

end

function Base.show(io::IO, x::MalmquistDEAModel)
    compact = get(io, :compact, false)

    Prod = prodchange(x, :Prod)
    EC   = prodchange(x, :EC)
    TC   = prodchange(x, :TC)
    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    periods = nperiods(x)

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
            show(io, CoefTable(hcat(Prod, EC, TC), ["M", "EC", "TC"], ["$i" for i in 1:n]))
        end
        if periods > 2
            show(io, CoefTable(hcat(Prod, EC, TC),
                               vcat(["M$i" for i in 1:2], ["EC$i" for i in 1:2], ["TC$i" for i in 1:2 ]),
                               ["$i" for i in 1:n]))
        end
        print(io, "\n")
        print(io, "M  = Malmquist Productivity Index \n")
        print(io, "EC = Efficiency Change \n")
        print(io, "TC = Technological Change")
    else

    end
end
