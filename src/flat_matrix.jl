# array of arrays but with a given number of arrays (M) and given total length (N)
# for internal use only
struct FlatMatrix{T,N,M}
    values::NTuple{N,T}
    indices::NTuple{M,Int64}

    function FlatMatrix(v::Vector{Vector{T}}) where {T}
        M = length(v)
        N = sum(length.(v))

        values = NTuple{N,T}(vcat(v...))
        indices = ntuple(i -> sum(length.(v[1:(i - 1)])), M)

        return new{T,N,M}(values, indices)
    end

    function FlatMatrix{T,N,M}(v::Vector{Vector{T}}) where {N,M,T}
        values = NTuple{N,T}(vcat(v...))
        indices = ntuple(i -> sum(length.(v[1:(i - 1)])), M)

        return new{T,N,M}(values, indices)
    end
end

function Base.getindex(m::FlatMatrix{T,N,M}, x::Int, y::Int) where {T,N,M}
    x <= M || throw(InvalidInputError("invalid indices ($x, $y) for flat matrix $m"))
    (x > 0 && y > 0) ||
        throw(InvalidInputError("invalid indices ($x, $y) for flat matrix $m"))
    if x == M
        m.indices[x] + y <= N
    else
        m.indices[x] + y <= m.indices[x + 1] ||
            throw(InvalidInputError("invalid indices ($x, $y) for flat matrix $m"))
    end
    return m.values[m.indices[x] + y]
end

function Base.length(m::FlatMatrix{T,N,M}, x::Int) where {T,N,M}
    (x <= M && x > 0) || throw(InvalidInputError("invalid index $x for flat matrix $m"))
    return x == M ? N - m.indices[x] : m.indices[x + 1] - m.indices[x]
end
