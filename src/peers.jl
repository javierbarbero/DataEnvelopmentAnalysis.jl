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
    dmuname::String
    dmunamesref::Vector{String}
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
    dmunames::Vector{String}
    dmunamesref::Vector{String}

    function DEAPeers(x::AbstractDEAModel; atol::Float64 = 1e-10, namesref::Vector{String} = Array{String}(undef, 0))
        if ! isdefined(x, :lambda)
            error("Model does not have info on peers.")
        end

        n = nobs(x)
        lambda = x.lambda
        nref = size(lambda, 2)
        dmunames = names(x)
        dmunamesref = String[]

        if n != size(lambda,1)
            error("Number of observation in the model and number of rows in peers matrix do not match")
        end

        # Remove very small numbers
        droptol!(lambda, atol)

        # If peers matrix is square
        if length(namesref) == 0
            if n == nref
                dmunamesref = dmunames
            else
                dmunamesref = ["Ref$i" for i in 1:nref]
            end
        elseif length(namesref) != nref
            error("Length of references names different to number of references DMUs")
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

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaio = dea(X, Y);

julia> peers(deaio)
DEA Peers
1: 1 ( 1.0 )
2: 4 ( 0.424978317432784 ) 7 ( 0.10928013876843023 )
3: 1 ( 1.134321653189578 ) 4 ( 0.43800539083557943 )
4: 4 ( 1.0 )
5: 4 ( 0.25738077214231636 ) 7 ( 0.04844814534443607 )
6: 7 ( 0.3333333333333333 )
7: 7 ( 1.0 )
8: 4 ( 1.0348650979425895 ) 7 ( 0.11457435012935832 )
9: 7 ( 1.1481481481481481 )
10: 4 ( 0.49056603773584906 ) 7 ( 0.4905660377358491 )
11: 11 ( 1.0 )
```
"""
peers(model::AbstractDEAModel; atol::Float64 = 1e-10, namesref::Vector{String} = Array{String}(undef, 0)) = DEAPeers(model, atol = atol, namesref = namesref) ;


"""
    ispeer(P::DEAPeers, i::Int64, j::Int64)
Return `true` if `j` is peer of decision making unit `i`.

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> P = peers(dea(X, Y));

julia> ispeer(P, 2, 4)
true
```
"""
function ispeer(P::DEAPeers, i::Int64, j::Int64)
    return P.lambda[i, j] != 0
end

"""
    ispeer(P::DEAPeers, i::String, j::String)
Return `true` if `j` is peer of decision making unit `i`.

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> firms = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; "I"; "J"; "K"]

julia> P = peers(dea(X, Y, names = firms));

julia> ispeer(P, "B", "D")
true
```
"""
function ispeer(P::DEAPeers, i::String, j::String)
    # Get corresponding id's of the names
    ival = findall(x -> x == i, P.dmunames)
    jval = findall(x -> x == j, P.dmunames)

    if length(ival) == 0
        error("Name $i does not exists")
    end
    if length(jval) == 0
        error("Name $j does not exists")
    end

    if length(ival) > 1
        error("Name $i is repeated. Search by name not possible.")
    end
    if length(jval) > 1
        error("Name $j is repeated. Search by name not possible.")
    end

    ival = ival[1]
    jval = jval[1]

    return ispeer(P, ival, jval)
end

"""
    ispeer(P::DEAPeersDMU, j::Int64)
Return `true` if `j` is peer.

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> P = peers(dea(X, Y));

julia> P2 = P[2];

julia> ispeer(P2, 4)
true
```
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

# Examples
```jldoctest
julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> firms = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; "I"; "J"; "K"]

julia> P = peers(dea(X, Y, names = firms));

julia> P2 = P[2]

julia> ispeer(P2, "D")
true
```
"""
function ispeer(p::DEAPeersDMU, j::String)::Bool
    return any(p.dmunamesref .== j)
end

# CONVERT FUNCTIONS

Base.convert(::Type{Matrix}, x::DEAPeers) = convert(Matrix, x.lambda);

Base.convert(::Type{SparseMatrixCSC}, x::DEAPeers) = x.lambda;
