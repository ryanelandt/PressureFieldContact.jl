
findMinSVSV(a::SVector{3,Float64}, b::SVector{3,Float64}) = min.(a, b)
findMinSVSV(sv::SVector{3, SVector{3,Float64}}) = min.(min.(sv[1], sv[2]),      sv[3])
findMinSVSV(sv::SVector{4, SVector{3,Float64}}) = min.(min.(sv[1], sv[2]), min.(sv[3], sv[4]))

findMaxSVSV(a::SVector{3,Float64}, b::SVector{3,Float64}) = max.(a, b)
findMaxSVSV(sv::SVector{3, SVector{3,Float64}}) = max.(max.(sv[1], sv[2]),      sv[3])
findMaxSVSV(sv::SVector{4, SVector{3,Float64}}) = max.(max.(sv[1], sv[2]), max.(sv[3], sv[4]))

function minMaxToCenterExtent(box_min::SVector{3,Float64}, box_max::SVector{3,Float64})
    center = (box_max + box_min) * 0.5
    extent = (box_max - box_min) * 0.5
    return center, extent
end

### calc_min_max
# calc_min_max(a::AABB) = a.c - a.e, a.c + a.e
function calc_min_max(a::OBB)
    δ = abs.(a.R) * a.e
    return a.c - δ, a.c + δ
end
function calc_min_max(point::Vector{SVector{3,T}}) where {T}
    min_val = SVector{3,Float64}(+Inf, +Inf, +Inf)
    max_val = SVector{3,Float64}(-Inf, -Inf, -Inf)
    for k = 1:length(point)
        point_k = point[k]
        min_val = min.(min_val, point_k)
        max_val = max.(max_val, point_k)
    end
    return min_val, max_val
end
function calc_min_max(a::SVector{3,Float64}, b::SVector{3,Float64})
    min_ = findMinSVSV(a, b)
    max_ = findMaxSVSV(a, b)
    return min_, max_
end

# ### calc_aabb
# function calc_aabb(vert::SVector{N, SVector{3,Float64}}) where {N}
#     box_min = findMinSVSV(vert)
#     box_max = findMaxSVSV(vert)
#     return calc_aabb(box_min, box_max)
# end
# function calc_aabb(arg_in)
#     min_val, max_val = calc_min_max(arg_in)
#     return calc_aabb(min_val, max_val)
# end
# function calc_aabb(a::SVector{3,Float64}, b::SVector{3,Float64})
#     min_val, max_val = calc_min_max(a, b)
#     center, extent = minMaxToCenterExtent(min_val, max_val)
#     return AABB(center, extent)
# end

### calc_obb
function calc_obb(vert::SVector{N, SVector{3,Float64}}) where {N}
    box_min = findMinSVSV(vert)
    box_max = findMaxSVSV(vert)
    return calc_obb(box_min, box_max)
end
function calc_obb(arg_in)
    min_val, max_val = calc_min_max(arg_in)
    return calc_obb(min_val, max_val)
end
function calc_obb(a::SVector{3,Float64}, b::SVector{3,Float64})
    min_val, max_val = calc_min_max(a, b)
    center, extent = minMaxToCenterExtent(min_val, max_val)
    return OBB(center, extent, one(SMatrix{3,3,Float64,9}))
end


function sortEdgeFace(v::SVector{3,Int64}, k::Int64)
    three = 3
    (1 <= k <= three) || error("a triangle has three sides")
    i1 = v[mod1(k + 1, three)]
    i2 = v[mod1(k + 2, three)]
    return SVector{2,Int64}(minmax(i1, i2))
end

function sortEdgeFace(v::SVector{4,Int64}, k::Int64)
    four = 4
    (1 <= k <= four) || error("a tetrahedron has four sides")
    i1 = v[mod1(k + 1, four)]
    i2 = v[mod1(k + 2, four)]
    i3 = v[mod1(k + 3, four)]
    i1, i2 = minmax(i1, i2)  # --> i2 is not lowest
    i2, i3 = minmax(i2, i3)  # --> i3 is highest
    i1, i2 = minmax(i1, i2)  # --> i1 and i2 are sorted
    return SVector{3,Int64}(i1, i2, i3)
end

function sort_so_big_ϵ_last(ϵ::SVector{4,Float64}, thing::SVector{4,T}) where {T}
	ϵ = abs.(ϵ)
	v, i = findmax(ϵ)
	return thing[tet_perm_by_num(i)]
end
