# This file contains types and structures for DEA peers
"""
    AbstractDEAPeersDMU
An abstract type representing the DEA Peers of a DMU.
"""
abstract type AbstractDEAPeersDMU end

"""
    DEAPeersDMU
An data structure representing the DEA Peers of a DMU.
"""
struct DEAPeersDMU <: AbstractDEAPeersDMU
    i::Int64
    J::Vector{Int64}
    V::Vector{Float64}
    dmuname::AbstractString
    dmunamesref::Vector{AbstractString}
end


Base.size(P::DEAPeersDMU) = (length(P.J),)

Base.length(P::DEAPeersDMU) = length(P.J)

Base.eltype(P::DEAPeersDMU) = Tuple{Tuple{Int64,String},Float64}

Base.firstindex(P::DEAPeersDMU) = 1

Base.lastindex(P::DEAPeersDMU) = length(P)

Base.getindex(P::DEAPeersDMU, i::Int) = return (P.J[i], P.dmunamesref[i]), P.V[i]

Base.iterate(P::DEAPeersDMU, state=1) = state > length(P) ? nothing : (P[state], state + 1)

function Base.show(io::IO, x::DEAPeersDMU)
    compact = get(io, :compact, false)

    J = x.J
    V = x.V
    dmuname = x.dmuname
    dmunamesref = x.dmunamesref

    print(io, dmuname, ": ")

    for p in x
        print(io, p[1][2], " ( ", p[2], " ) ")
    end

end


"""
    AbstractDEAPeers
An abstract type representing a DEA Peers.
"""
abstract type AbstractDEAPeers end

"""
    DEAPeers
An data structure representing a DEA Peers.
"""
struct DEAPeers <: AbstractDEAPeers
    n::Int64
    nref::Int64
    lambda::SparseMatrixCSC{Float64, Int64}
    I::Vector{Int64}
    J::Vector{Int64}
    V::Vector{Float64}
    dmunames::Vector{AbstractString}
    dmunamesref::Vector{AbstractString}

    function DEAPeers(x::AbstractDEAModel; atol::Float64 = 1e-10, namesref::Union{Vector{<: AbstractString},Nothing} = nothing)
        if ! isdefined(x, :lambda)
            throw(ArgumentError("Model does not have info on peers"));
        end

        n = nobs(x)
        lambda = deepcopy(x.lambda)
        nref = size(lambda, 2)
        dmunames = names(x)
        dmunamesref = String[]

        # Remove very small numbers
        droptol!(lambda, atol)

        # If peers matrix is square
        if namesref === nothing
            if n == nref
                dmunamesref = dmunames
            else
                dmunamesref = ["Ref$i" for i in 1:nref]
            end
        elseif length(namesref) != nref
            throw(DimensionMismatch("length of `namesref` and number of reference DMUs are not equal"));
        else
            dmunamesref = namesref
        end

        # Extract I, J and V vectors from SparseArray and order by I and J
        I, J, V = findnz(lambda)
        npeers = length(I)
        IJV = [I J V]
        IJV = sortslices(IJV, dims = 1, by = x -> (x[1], x[2]))

        I = convert.(Int, IJV[:,1])
        J = convert.(Int, IJV[:,2])
        V = IJV[:,3]

        new(n, nref, lambda, I, J, V, dmunames, dmunamesref)

    end

end

Base.size(P::DEAPeers) = (P.n,)

Base.length(P::DEAPeers) = P.n

Base.eltype(P::DEAPeers) = DEAPeersDMU

Base.firstindex(P::DEAPeers) = 1

Base.lastindex(P::DEAPeers) = length(P)

function Base.getindex(P::DEAPeers, i::Int)::DEAPeersDMU

    I = P.I
    J = P.J
    V = P.V
    dmunames = P.dmunames
    dmunamesref = P.dmunamesref

    DEAPeersDMU(i, J[I .== i], V[I .== i], dmunames[i], dmunamesref[unique(J[I .== i])])

end

Base.iterate(P::DEAPeers, state = 1) = state > length(P) ? nothing : (P[state], state+1)

function Base.show(io::IO, x::DEAPeers)
    compact = get(io, :compact, false)

    n = length(x)

    println(io, "DEA Peers")

    for p in x
        println(io, p)
    end

end

"""
    peers(model::AbstractDEAModel)
Return peers of a DEA model.

# Optional Arguments
- `atol=1e-10`: tolerance for zero values.
- `namesref`: a vector of strings with the names of the decision making units in the reference set.
"""
peers(model::AbstractDEAModel; atol::Float64 = 1e-10, namesref::Union{Vector{<: AbstractString},Nothing} = nothing) = DEAPeers(model, atol = atol, namesref = namesref) ;

"""
    peersmatrix(model::AbstractDEAModel)
Return peers matrix of a DEA model.
"""
peersmatrix(model::AbstractDEAModel) = model.lambda ;

"""
    ispeer(P::DEAPeers, i::Int64, j::Int64)
Return `true` if `j` is peer of decision making unit `i`.
"""
function ispeer(P::DEAPeers, i::Int64, j::Int64)
    return P.lambda[i, j] != 0
end

"""
    ispeer(P::DEAPeers, i::String, j::String)
Return `true` if `j` is peer of decision making unit `i`.
"""
function ispeer(P::DEAPeers, i::String, j::String)
    # Get corresponding id's of the names
    ival = findall(x -> x == i, P.dmunames)
    jval = findall(x -> x == j, P.dmunames)

    if length(ival) == 0
        throw(ArgumentError("Name $i does not exists"));
    end
    if length(jval) == 0
        throw(ArgumentError("Name $j does not exists"));
    end

    if length(ival) > 1
        throw(ArgumentError("Name $i is repeated. Search by name not possible."));
    end
    if length(jval) > 1
        throw(ArgumentError("Name $j is repeated. Search by name not possible."));
    end

    ival = ival[1]
    jval = jval[1]

    return ispeer(P, ival, jval)
end

"""
    ispeer(P::DEAPeersDMU, j::Int64)
Return `true` if `j` is peer.
"""
function ispeer(p::DEAPeersDMU, j::Int64)::Bool
    if j <= 0
        throw(BoundsError(p, j))
    end

    return any(p.J .== j)
end

"""
    ispeer(P::DEAPeers, j::String)
Return `true` if `j` is peer.
"""
function ispeer(p::DEAPeersDMU, j::String)::Bool
    return any(p.dmunamesref .== j)
end

"""
    sum(P::DEAPeers; dims)

Sum elements of the peers matrix over the given dimension.
"""
Base.sum(P::DEAPeers; dims) = sum(P.lambda; dims = dims)

# CONVERT FUNCTIONS

Base.convert(::Type{Matrix}, x::DEAPeers) = convert(Matrix, x.lambda);

Base.convert(::Type{SparseMatrixCSC}, x::DEAPeers) = x.lambda;
