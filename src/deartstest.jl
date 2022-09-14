
# This file contains functions for the DEA Returns to Scale Test
struct DEAReturnsToScaleTest
    nreps::Int        
    S::Float64
    SB::Vector{Float64}    
    p::Float64 
    h::Float64
end

"""
    deartstest(X, Y)
Compute the DEA Returns to Scale (RTS) test using the bootstrap radial model for inputs X and outputs Y.

# Optional Arguments
- `nreps=200`: number of bootstrap replications.
- `rng=default_rng()`: random number generator.
- `orient=:Input`: chooses the radially oriented input mode. For the radially oriented output model choose `:Output`.
"""
function deartstest(X::Union{Matrix,Vector}, Y::Union{Matrix,Vector}; 
    orient::Symbol = :Input, nreps::Int = 200, rng::AbstractRNG = default_rng())::DEAReturnsToScaleTest

    # Observed S
    crs = efficiency(dea(X, Y, rts = :CRS, orient = orient))
    vrs = efficiency(dea(X, Y, rts = :VRS, orient = orient))

    S = 0
    if orient == :Input
        S = sum(crs) / sum(vrs)
    else
        S = sum(1 ./ crs) / sum(1 ./ vrs)
    end

    # Bootstrapped S
    rng2 = deepcopy(rng)

    crsboot = deaboot(X, Y, orient = orient, rts = :CRS, nreps = nreps, rng = rng)
    crsB = crsboot.effB

    vrsboot = deaboot(X, Y, orient = orient, rts = :VRS, nreps = nreps, rng = rng2, effref = crs)
    vrsB = vrsboot.effB

    if orient == :Input
        SB = sum(crsB, dims = 1) ./ sum(vrsB, dims = 1)
    else
        SB = sum(1 ./ crsB, dims = 1) ./ sum(1 ./ vrsB, dims = 1)
    end
    SB = vec(SB)

    lower = sum(SB .< S)
    p = (lower + 1) / nreps

    return DEAReturnsToScaleTest(nreps, S, SB, p, crsboot.h)
end

"""
    criticalvalue(model::BootstrapRadialDEAModel, alpha = 0.05)
Return the critical value of the DEA Returns to Scale Test, with level `alpha`.
"""
function criticalvalue(x::DEAReturnsToScaleTest, alpha ::Float64 = 0.05)
    SB = x.SB
    nreps = x.nreps

    SBsorted = sort(SB)
    critval = mean( SBsorted[Int.([floor(alpha * nreps), 
                                   ceil(alpha * nreps)])] )
    return critval
end

function Base.show(io::IO, x::DEAReturnsToScaleTest)
    compact = get(io, :compact, false)

    nreps = x.nreps    
    S = x.S
    p = x.p
    h = x.h
    alpha = 0.05
    critval = criticalvalue(x, alpha)

    if !compact
        print(io, "DEA Returns to Scale (RTS) Test \n")
        print(io, "--------------------------------\n")
        print(io, "\n")
        print(io, "  H0: Globally CRS \n")
        print(io, "  H1: VRS \n")
        print(io, "\n")
        print(io, "  Bootstrap replications: ", nreps, "\n")
        print(io, "  Bandwidth = ", round(h, digits = 5), "\n")
        print(io, "\n")
        print(io, "  Scale efficiency: ", round(S, digits = 4), "\n")
        print(io, "  Critical value (α = $alpha): ", round(critval, digits = 4), "\n")
        print(io, "  p-value: ", round(p, digits = 4), "\n")
    end

end
