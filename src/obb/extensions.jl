
NumericalTricks.area(a::BB)   where {BB <: BoundingBox} = 8 * dot(a.e, SVector{3,Float64}(a.e[2], a.e[3], a.e[1]))
NumericalTricks.volume(a::BB) where {BB <: BoundingBox} = 8 * prod(a.e)
