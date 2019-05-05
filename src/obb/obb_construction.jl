
function tet_perm_by_num(n::Int64)
    (n == 1) && (return SVector{4,Int64}(2,4,3,1))
    (n == 2) && (return SVector{4,Int64}(4,1,3,2))
    (n == 3) && (return SVector{4,Int64}(1,4,2,3))
    (n == 4) && (return SVector{4,Int64}(1,2,3,4))
    error("n must be 1:4")
end

get_projections(p::SVector{3,SVector{3,Float64}}, ê) = SVector(dot(p[1], ê), dot(p[2], ê), dot(p[3], ê))
get_projections(p::SVector{4,SVector{3,Float64}}, ê) = SVector(dot(p[1], ê), dot(p[2], ê), dot(p[3], ê), dot(p[4], ê))

function make_obb(p::SVector{N,SVector{3,Float64}}, i_start::Int64) where {N}
    (N == 3) || (N == 4) || error("must have 3 or 4 points")
    ê_1 = normalize(p[mod1(i_start + 1, 3)] - p[i_start])
    ê_3 = triangleNormal(p[1], p[2], p[3])
    ê_2 = cross(ê_3, ê_1)
    proj_1 = get_projections(p, ê_1)
    proj_2 = get_projections(p, ê_2)
    proj_3 = get_projections(p, ê_3)
    p_min = SVector(minimum(proj_1), minimum(proj_2), minimum(proj_3))
    p_max = SVector(maximum(proj_1), maximum(proj_2), maximum(proj_3))
    c, e = minMaxToCenterExtent(p_min, p_max)
    R = hcat(ê_1, ê_2, ê_3)
    return OBB(R * c, e, R)
end

fit_tri_obb(p::SVector{3,SVector{3,Float64}}) = make_obb(p, 1)  # all bounding boxes have the same surface area
function fit_tet_obb(p::SVector{4,SVector{3,Float64}}, ϵ_tet::SVector{4,Float64})
    (0.0 < volume(p)) || error("inverted tet")
    p = sort_so_big_ϵ_last(ϵ_tet, p)
    obb_1 = make_obb(p, 1)
    obb_2 = make_obb(p, 2)
    obb_3 = make_obb(p, 3)
    area_1 = area(obb_1)
    area_2 = area(obb_2)
    area_3 = area(obb_3)
    (max(        area_2, area_3) <= area_1) && (return obb_1)
    (max(area_1,         area_3) <= area_2) && (return obb_2)
    (max(area_1, area_2        ) <= area_3) && (return obb_3)
end
