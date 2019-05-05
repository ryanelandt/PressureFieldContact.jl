
abstract type BoundingBox end

struct OBB <: BoundingBox
    c::SVector{3,Float64}
    e::SVector{3,Float64}
    R::SMatrix{3,3,Float64,9}
    OBB(c::SVector{3,Float64}, e::SVector{3,Float64}, R::SMatrix{3,3,Float64,9}) = new(c, e, R)
end

function (::Type{BB_Type})(a::BB_Type, b::BB_Type) where {BB_Type <: BoundingBox}
    min_1, max_1 = calc_min_max(a)
    min_2, max_2 = calc_min_max(b)
    return calc_obb(SVector{4,SVector{3,Float64}}(min_1, max_1, min_2, max_2))
end
