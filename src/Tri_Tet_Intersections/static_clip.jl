
"""
Calls clip
"""
function clip_in_tet_coordinates(p::poly_eight{4,T}) where {T}
    if length(p) == 3
        return clip_in_tet_coordinates(p.v[1], p.v[2], p.v[3])
    elseif length(p) == 4
        return clip_in_tet_coordinates(p.v[1], p.v[2], p.v[3], p.v[4])
    else
        error("something is wrong")
    end
end

function clip_in_tet_coordinates(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}) where {T}
    return clip(z1, z2, z3, 1)
end

function clip_in_tet_coordinates(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}) where {T}
    return clip(z1, z2, z3, z4, 1)
end

"""
Together with cut_clip, implements the clipping of polygons by tetrahedra using the Sutherland–Hodgman algorithm.
The Sutherland-Hodgman algorithm essentially clips a polygon by one plane at a time.
You should not call this function directly.
Clipped polygons are assumed to have at most 8 sides.
This implementation should correctly handle the case where 2 of fewer vertices lie exactly on a plane.
"""
function clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, i::Int64) where {T}
    # Polygons are represented as vertices in tetrahedral coordinates.
    # If all components of a vertex are positive this means that the vertex lies within the tetrahedron
    # Negative components mean that a vertix lies outside a face.
    # This function permutes the vertices until the first is non-positive and the second is positive for cut_clip.

    (i == 5) && (return poly_eight(3, (z1, z2, z3, z1, z1, z1, z1, z1)))  # there is no 5th plane return polygon

    s = SVector{3,T}(z1[i], z2[i], z3[i])
    is_non_pos = (s .<= 0.0)
    all(is_non_pos) && (return poly_eight{4,T}())  # entire polygon clipped away
    if all(0.0 .<= s)  # non_negative
        return clip(z1, z2, z3, i+1)  # all points inside this plane go to next plane
    else
        is_non_pos[1] && !is_non_pos[2] && (return cut_clip(z1, z2, z3, i))
        is_non_pos[2] && !is_non_pos[3] && (return cut_clip(z2, z3, z1, i))
        is_non_pos[3] && !is_non_pos[1] && (return cut_clip(z3, z1, z2, i))
    end
    error("Non-finite vertex likely")
end

function clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, i::Int64) where {T}
    (i == 5) && (return poly_eight(4, (z1, z2, z3, z4, z1, z1, z1, z1)))

    s = SVector{4,T}(z1[i], z2[i], z3[i], z4[i])
    is_non_pos = (s .<= 0.0)
    all(is_non_pos) && (return poly_eight{4,T}())
    if all(0.0 .<= s)  # non_negative
        return clip(z1, z2, z3, z4, i + 1)
    else
        is_non_pos[1] && !is_non_pos[2] && (return cut_clip(z1, z2, z3, z4, i))
        is_non_pos[2] && !is_non_pos[3] && (return cut_clip(z2, z3, z4, z1, i))
        is_non_pos[3] && !is_non_pos[4] && (return cut_clip(z3, z4, z1, z2, i))
        is_non_pos[4] && !is_non_pos[1] && (return cut_clip(z4, z1, z2, z3, i))
    end
    error("Non-finite vertex likely")
end

function clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, z5::SVector{4,T}, i::Int64) where {T}
    (i == 5) && (return poly_eight(5, (z1, z2, z3, z4, z5, z1, z1, z1)))

    s = SVector{5,T}(z1[i], z2[i], z3[i], z4[i], z5[i])
    is_non_pos = (s .<= 0.0)
    all(is_non_pos) && (return poly_eight{4,T}())
    if all(0.0 .<= s)  # non_negative
        return clip(z1, z2, z3, z4, z5, i + 1)
    else
        is_non_pos[1] && !is_non_pos[2] && (return cut_clip(z1, z2, z3, z4, z5, i))
        is_non_pos[2] && !is_non_pos[3] && (return cut_clip(z2, z3, z4, z5, z1, i))
        is_non_pos[3] && !is_non_pos[4] && (return cut_clip(z3, z4, z5, z1, z2, i))
        is_non_pos[4] && !is_non_pos[5] && (return cut_clip(z4, z5, z1, z2, z3, i))
        is_non_pos[5] && !is_non_pos[1] && (return cut_clip(z5, z1, z2, z3, z4, i))
    end
    error("Non-finite vertex likely")
end

function clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, z5::SVector{4,T}, z6::SVector{4,T}, i::Int64) where {T}
    (i == 5) && (return poly_eight(6, (z1, z2, z3, z4, z5, z6, z1, z1)))


    s = SVector{6,T}(z1[i], z2[i], z3[i], z4[i], z5[i], z6[i])
    is_non_pos = (s .<= 0.0)
    all(is_non_pos) && (return poly_eight{4,T}())
    if all(0.0 .<= s)  # non_negative
        return clip(z1, z2, z3, z4, z5, z6, i + 1)
    else
        is_non_pos[1] && !is_non_pos[2] && (return cut_clip(z1, z2, z3, z4, z5, z6, i))
        is_non_pos[2] && !is_non_pos[3] && (return cut_clip(z2, z3, z4, z5, z6, z1, i))
        is_non_pos[3] && !is_non_pos[4] && (return cut_clip(z3, z4, z5, z6, z1, z2, i))
        is_non_pos[4] && !is_non_pos[5] && (return cut_clip(z4, z5, z6, z1, z2, z3, i))
        is_non_pos[5] && !is_non_pos[6] && (return cut_clip(z5, z6, z1, z2, z3, z4, i))
        is_non_pos[6] && !is_non_pos[1] && (return cut_clip(z6, z1, z2, z3, z4, z5, i))
    end
    error("Non-finite vertex likely")
end

function clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, z5::SVector{4,T}, z6::SVector{4,T}, z7::SVector{4,T}, i::Int64) where {T}
    (i == 5) && (return poly_eight(7, (z1, z2, z3, z4, z5, z6, z7, z1)))

    s = SVector{7,T}(z1[i], z2[i], z3[i], z4[i], z5[i], z6[i], z7[i])
    is_non_pos = (s .<= 0.0)
    all(is_non_pos) && (return poly_eight{4,T}())
    if all(0.0 .<= s)  # non_negative
        return clip(z1, z2, z3, z4, z5, z6, z7, i + 1)
    else
        is_non_pos[1] && !is_non_pos[2] && (return cut_clip(z1, z2, z3, z4, z5, z6, z7, i))
        is_non_pos[2] && !is_non_pos[3] && (return cut_clip(z2, z3, z4, z5, z6, z7, z1, i))
        is_non_pos[3] && !is_non_pos[4] && (return cut_clip(z3, z4, z5, z6, z7, z1, z2, i))
        is_non_pos[4] && !is_non_pos[5] && (return cut_clip(z4, z5, z6, z7, z1, z2, z3, i))
        is_non_pos[5] && !is_non_pos[6] && (return cut_clip(z5, z6, z7, z1, z2, z3, z4, i))
        is_non_pos[6] && !is_non_pos[7] && (return cut_clip(z6, z7, z1, z2, z3, z4, z5, i))
        is_non_pos[7] && !is_non_pos[1] && (return cut_clip(z7, z1, z2, z3, z4, z5, z6, i))
    end
    error("Non-finite vertex likely")
end

"""
Removes a non-positive vertex from a polygon.
"""
function cut_clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, i::Int64) where {T}
    # z1 is gauranteed to be non-positive
    # z2 is gauranteed to be positive

    z_start = clip_node(z1, z2, i)
    if 0.0 < z3[i]
        z_end = clip_node(z1, z3, i)
        return clip(z_start, z2, z3, z_end, i + 1)
    else
        z_end = clip_node(z3, z2, i)
        return clip(z_start, z2, z_end, i + 1)
    end
end

function cut_clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, i::Int64) where {T}
    (z3[i] <= 0.0) && (return cut_clip(z1, z2, z3, i))
    z_start = clip_node(z1, z2, i)
    if 0.0 < z4[i]
        z_end = clip_node(z1, z4, i)
        return clip(z_start, z2, z3, z4, z_end, i + 1)
    else
        z_end = clip_node(z4, z3, i)
        return clip(z_start, z2, z3, z_end, i + 1)
    end
end

function cut_clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, z5::SVector{4,T}, i::Int64) where {T}
    (z4[i] <= 0.0) && (return cut_clip(z1, z2, z3, z4, i))
    z_start = clip_node(z1, z2, i)
    if 0.0 < z5[i]
        z_end = clip_node(z1, z5, i)
        return clip(z_start, z2, z3, z4, z5, z_end, i + 1)
    else
        z_end = clip_node(z5, z4, i)
        return clip(z_start, z2, z3, z4, z_end, i + 1)
    end
end

function cut_clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, z5::SVector{4,T}, z6::SVector{4,T}, i::Int64) where {T}
    (z5[i] <= 0.0) && (return cut_clip(z1, z2, z3, z4, z5, i))
    z_start = clip_node(z1, z2, i)
    if 0.0 <= z6[i]
        z_end = clip_node(z1, z6, i)
        return clip(z_start, z2, z3, z4, z5, z6, z_end, i + 1)
    else
        z_end = clip_node(z6, z5, i)
        return clip(z_start, z2, z3, z4, z5, z_end, i + 1)
    end
end

function cut_clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, z5::SVector{4,T}, z6::SVector{4,T}, z7::SVector{4,T}, i::Int64) where {T}
    (z6[i] <= 0.0) && (return cut_clip(z1, z2, z3, z4, z5, z6, i))
    z_start = clip_node(z1, z2, i)
    if 0.0 <= z7[i]
        z_end = clip_node(z1, z7, i)
        return poly_eight(8, (z_start, z2, z3, z4, z5, z6, z7, z_end))
    else
        z_end = clip_node(z7, z6, i)
        return poly_eight(7, (z_start, z2, z3, z4, z5, z6, z_end, z_start))
    end
end

function clip_node(z_non::SVector{4,T}, z_pos::SVector{4,T}, i::Int64) where {T}
    v_non = z_non[i]
    v_pos = z_pos[i]
    return weightPoly(z_non, z_pos, v_non, v_pos)
end
