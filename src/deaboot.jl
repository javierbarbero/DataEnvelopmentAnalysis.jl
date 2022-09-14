# This file contains functions for the Bootstrap Radial DEA model
"""
AbstractBootstrapDEAModel
An abstract type representing a bootstrap DEA model.
"""
abstract type AbstractBootstrapDEAModel <: AbstractTechnicalDEAModel end

"""
BootstrapRadialDEAModel
An data structure representing a bootstrap radial DEA model.
"""
struct BootstrapRadialDEAModel <: AbstractBootstrapDEAModel
    n::Int64
    m::Int64
    s::Int64
    orient::Symbol
    rts::Symbol
    disposX::Symbol
    disposY::Symbol
    dmunames::Union{Vector{AbstractString},Nothing}
    eff::Vector
    effref::Vector
    effbias::Vector
    effB::Matrix
    h::Float64
end

"""
    deaboot(X, Y)
Compute the bootstrap radial model using data envelopment analysis for inputs X and outputs Y.

# Optional Arguments
- `nreps=200`: number of bootstrap replications.
- `rng=default_rng()`: random number generator.
- `orient=:Input`: chooses the radially oriented input mode. For the radially oriented output model choose `:Output`.
- `rts=:CRS`: chooses constant returns to scale. For variable returns to scale choose `:VRS`.
- `Xref=X`: Identifies the reference set of inputs against which the units are evaluated.
- `Yref=Y`: Identifies the reference set of outputs against which the units are evaluated.
- `disposX=:Strong`: chooses strong disposability of inputs. For weak disposability choose `:Weak`.
- `disposY=:Strong`: chooses strong disposability of outputs. For weak disposability choose `:Weak`.
- `names`: a vector of strings with the names of the decision making units.
"""
function deaboot(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector};
    nreps::Int = 200, rng::AbstractRNG = default_rng(), effref::Union{Vector,Nothing} = nothing,
    orient::Symbol = :Input, rts::Symbol = :CRS,
    Xref::Union{Matrix,Vector,Nothing} = nothing, Yref::Union{Matrix, Vector,Nothing} = nothing,
    disposX::Symbol = :Strong, disposY::Symbol = :Strong,
    names::Union{Vector{<: AbstractString},Nothing} = nothing,
    optimizer::Union{DEAOptimizer,Nothing} = nothing, progress::Bool = false)::BootstrapRadialDEAModel

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

    if disposX != :Strong && disposX != :Weak
        throw(ArgumentError("`disposX` must be :Strong or :Weak"));
    end

    if disposY != :Strong && disposY != :Weak
        throw(ArgumentError("`disposY` must be :Strong or :Weak"));
    end

    # Default optimizer
    if optimizer === nothing 
        optimizer = DEAOptimizer(:LP)
    end

    # Radial DEA Model
    n = nx

    if isnothing(effref)
        effref = efficiency(dea(X, Y, orient = orient, rts = rts, slack = false, 
            Xref = Xref, Yref = Yref, disposX = disposX, disposY = disposY, optimizer = optimizer))
    end

    h = dea_bandwidth(effref, orient)

    if orient == :Input
        effi = 1 ./ deepcopy(effref)
    else
        effi = deepcopy(effref)
    end

    # Boostrap replications
    effB = SharedMatrix{Float64}(n, nreps)
    fill!(effB, 0.0)

    rng2 = deepcopy(rng)
    eff2m = [2 .- effi; effi]

    sampB = zeros(n, nreps)
    randB = zeros(n, nreps)
    for b in 1:nreps
        sampB[:, b] = sample(rng, eff2m, n, replace = true)
        randB[:, b] = rand(rng2, Normal(0, 1), n);
    end
    
    p = progressbarmeter(nreps, desc = "Computing Bootsrap Radial DEA Model...", progress = progress)

    @sync @distributed for b in 1:nreps
        # 1. Get a sample from the 2m set
        effn = sampB[:, b]

        # 2. Perturbate
        effn_star = effn + h .* randB[:,b]

        # 3. Refine and correct for the mean and variance of smoothed values
        effn_star = mean(effn) .+ (effn .- mean(effn)) ./ (sqrt( 1 + h^2 / var(eff2m)));

        # 4 .Reflect values
        effn_starr = effn_star
        effn_starr[effn_star .< 1] = 2 .- effn_star[effn_star .< 1]

        # Generate inefficient inputs or outputs and perform DEA
        if orient == :Input
            effn_starr = 1 ./ effn_starr
            Xrefb = repeat(effref ./ effn_starr, 1, m) .* Xref
            Yrefb = Yref
        else
            Xrefb = Xref
            Yrefb = repeat(effref ./ effn_starr, 1, s) .* Yref
        end

        effB[:, b] = efficiency(dea(X, Y, orient = orient, rts = rts, slack = false, 
            Xref = Xrefb, Yref = Yrefb, disposX = disposX, disposY = disposY, optimizer = optimizer))

        if orient == :Input
            effB[:, b] = 1 ./ effB[:, b]
        end

        next!(p)
    end

    # Bootstrap efficiency
    effbias = vec(mean(effB, dims = 2) .- effi)
    effbc = vec(effi .- effbias)

    # Invert if input oriented
    if orient == :Input
        effB = 1 ./ effB
        effbc = 1 ./ effbc
        effbias = effref .- effbc
    end

    return BootstrapRadialDEAModel(n, m, s, orient, rts, disposX, disposY, names, effbc, effref, effbias, effB, h)

end

function dea_bandwidth(x::Vector{Float64}, orient::Symbol)::Float64
    # Get only inefficient units
    effn = x[.! isapprox.(x, 1.0, atol = 1e-5)]

    n = length(x)
    m = length(effn)

    # Build the 2m set (reflected data)
    eff2m = [2 .- effn; effn]
    
    s2m = std(eff2m)
    r2m = iqr(eff2m)

    # Compute the optimmal bandwidth (After equation 3.24)
    hm = 0.9 * minimum([s2m, r2m ./ 1.349]) * (2 * m)^(-1/5);

    # Adjust bandwidth for scale and sample size (Equation 3.26) 
    if orient == :Input
        sx = std( 1 ./ x)
    else
        sx = std(x)
    end

    h = hm * ((2 * m) / n)^(1/5) * (sx / s2m ); 

    return h
end

"""
    confint(model::BootstrapRadialDEAModel; level::Real=0.95)
Compute confidence intervals for efficiency scores, with confidence level `level` (by default 95%).
"""
function confint(model::BootstrapRadialDEAModel; level::Real=0.95)
    n = nobs(model)
    orient = model.orient
    effref = model.effref
    effB = model.effB
    alpha = 1 - level

    effc = zeros(n, 2)
    for i in 1:n
        if orient == :Input
            effc[i, 1:2] = (1 ./ effref[i]) .+ quantile(( 1 ./ effref[i]) .- (1 ./ effB[i,:]), [0.5 * alpha, 1 - 0.5 * alpha])
        else
            effc[i, 1:2] = effref[i] .+ quantile(effref[i] .- effB[i,:], [0.5 * alpha, 1 - 0.5 * alpha])
        end
    end

    if orient == :Input
        effc = 1 ./ effc
        effc = effc[:, [2, 1]]
    end

    return effc
end

"""
    bandwidth(model::BootstrapRadialDEAModel)
Return the optimal bandwidth of a bootstrap DEA model.
"""
function bandwidth(model::BootstrapRadialDEAModel)
    return model.h
end

"""
    bias(model::BootstrapRadialDEAModel)
Return the bias from bootstrap DEA model.
"""
function bias(model::BootstrapRadialDEAModel)
    return model.effbias
end

function Base.show(io::IO, x::BootstrapRadialDEAModel)
    compact = get(io, :compact, false)

    n = nobs(x)
    m = ninputs(x)
    s = noutputs(x)
    disposX = x.disposX
    disposY = x.disposY
    dmunames = names(x)

    effref = x.effref
    eff = efficiency(x)
    effbias = bias(x)
    effc = confint(x, level = 0.95)
    effL = effc[:, 1]
    effU = effc[:, 2]

    if !compact
        print(io, "Bootstrap Radial DEA Model \n")
        print(io, "DMUs = ", n)
        print(io, "; Inputs = ", m)
        print(io, "; Outputs = ", s)
        print(io, "\n")
        print(io, "Orientation = ", string(x.orient))
        print(io, "; Returns to Scale = ", string(x.rts))
        print(io, "\n")
        print(io, "Bootstrap replications = ", size(x.effB, 2))
        print(io, "\n")
        if disposX == :Weak print(io, "Weak disposability of inputs \n") end
        if disposY == :Weak print(io, "Weak disposability of outputs \n") end

        show(io, CoefTable(hcat(effref, eff, effbias, effL, effU), ["Reference", "Corrected", "Bias", "Lower 95%", "Upper 95%"], dmunames))   
        
        print(io, "\nBandwidth = ", round(x.h, digits = 5))
    end

end
