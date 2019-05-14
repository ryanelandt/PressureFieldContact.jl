
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


function sort_so_big_ϵ_last(ϵ::SVector{4,Float64}, thing::SVector{4,T}) where {T}
	ϵ = abs.(ϵ)
	v, i = findmax(ϵ)
	return thing[tet_perm_by_num(i)]
end
