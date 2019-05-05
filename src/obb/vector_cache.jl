mutable struct VectorCache{T}
    ind_fill::Int64
    ind_max::Int64
    vec::Vector{T}
    function VectorCache{T}() where {T}
        isbitstype(T) || error("This cache requires an isbits type.")
        ind_max = 64
        new{T}(-9999, ind_max, Vector{T}(undef, ind_max))
    end
end

function expand!(vc::VectorCache{T}) where {T}  # TODO: make this function more elegant
    vc.ind_max += vc.ind_max
    resize!(vc.vec, vc.ind_max)
    return nothing
end

function returnNext(vc::VectorCache{T}) where {T}
    (vc.ind_fill == vc.ind_max) && expand!(vc)
    vc.ind_fill += 1
    return vc.vec[vc.ind_fill]
end

function addCacheItem!(vc::VectorCache{T}, item_T::T) where {T}
  vc.ind_fill += 1
  (vc.ind_max < vc.ind_fill) && expand!(vc)
  vc.vec[vc.ind_fill] = item_T
  return nothing
end

function Base.empty!(vc::VectorCache{T}) where {T}
  vc.ind_fill = 0
  return nothing
end

Base.isempty(vc::VectorCache{T}) where {T} = (vc.ind_fill == 0)
Base.length(vc::VectorCache{T}) where {T} = vc.ind_fill
Base.@propagate_inbounds Base.getindex(vc::VectorCache{T}, i::Int) where {T} = vc.vec[i]
