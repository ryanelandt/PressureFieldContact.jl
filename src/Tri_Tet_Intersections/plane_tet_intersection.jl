
"""
Clips a plane by a tetrahedron.
The plane takes the form of a static row matrix
The tet takes the form of a static matrix
"""
function clip_plane_tet(plane::SMatrix{1,4,T,4}, tet::SMatrix{4,4,T,16}) where {T}
    function tet_mat_to_tuple(tet::SMatrix{4,4,T,16}) where {T}
        v_1 = SVector{3,T}(tet[1],  tet[2],  tet[3])
        v_2 = SVector{3,T}(tet[5],  tet[6],  tet[7])
        v_3 = SVector{3,T}(tet[9],  tet[10], tet[11])
        v_4 = SVector{3,T}(tet[13], tet[14], tet[15])
        return (v_1, v_2, v_3, v_4)
    end

    proj = plane * tet

    bool_neg = proj .< 0.0
    n_neg = sum(bool_neg)

    bool_pos = 0.0 .< proj
    n_pos = sum(bool_pos)

    if (n_pos == 0) || (n_neg == 0)
        return poly_eight{3,T}()
    else
        v = tet_mat_to_tuple(tet)
        if n_pos == 1
            bool_pos[1] && (return clip_plane_tet_1(proj, v))
            bool_pos[2] && (return clip_plane_tet_2(proj, v))
            bool_pos[3] && (return clip_plane_tet_3(proj, v))
            bool_pos[4] && (return clip_plane_tet_4(proj, v))
        elseif n_neg == 1
            bool_neg[1] && (return clip_plane_tet_1(proj, v))
            bool_neg[2] && (return clip_plane_tet_2(proj, v))
            bool_neg[3] && (return clip_plane_tet_3(proj, v))
            bool_neg[4] && (return clip_plane_tet_4(proj, v))
        else  # if (n_pos == 2) && (n_neg == 2)
            (bool_pos[1] == bool_pos[2]) && (return clip_plane_tet_12(proj, v))
            (bool_pos[1] == bool_pos[3]) && (return clip_plane_tet_13(proj, v))
            (bool_pos[1] == bool_pos[4]) && (return clip_plane_tet_14(proj, v))
        # else
        #     error("something is wrong")
        end
    end
end

function proj_weight_poly(proj::SMatrix{1,4,T,4}, v::NTuple{4,SVector{3,T}}, i_1::Int64, i_2::Int64) where {T}
    return weightPoly(v[i_1], v[i_2], proj[i_1], proj[i_2])
end

function clip_plane_tet_1(proj::SMatrix{1,4,T,4}, v::NTuple{4,SVector{3,T}}) where {T}
    a = proj_weight_poly(proj, v, 2, 1)
    b = proj_weight_poly(proj, v, 4, 1)
    c = proj_weight_poly(proj, v, 3, 1)
    tup = 0.0 < proj[1] ? (a, b, c) : (c, b, a)
    return poly_eight(tup)
end
function clip_plane_tet_2(proj::SMatrix{1,4,T,4}, v::NTuple{4,SVector{3,T}}) where {T}
    a = proj_weight_poly(proj, v, 1, 2)
    b = proj_weight_poly(proj, v, 3, 2)
    c = proj_weight_poly(proj, v, 4, 2)
    tup = 0.0 < proj[2] ? (a, b, c) : (c, b, a)
    return poly_eight(tup)
end
function clip_plane_tet_3(proj::SMatrix{1,4,T,4}, v::NTuple{4,SVector{3,T}}) where {T}
    a = proj_weight_poly(proj, v, 1, 3)
    b = proj_weight_poly(proj, v, 4, 3)
    c = proj_weight_poly(proj, v, 2, 3)
    tup = 0.0 < proj[3] ? (a, b, c) : (c, b, a)
    return poly_eight(tup)
end
function clip_plane_tet_4(proj::SMatrix{1,4,T,4}, v::NTuple{4,SVector{3,T}}) where {T}
    a = proj_weight_poly(proj, v, 1, 4)
    b = proj_weight_poly(proj, v, 2, 4)
    c = proj_weight_poly(proj, v, 3, 4)
    tup = 0.0 < proj[4] ? (a, b, c) : (c, b, a)
    return poly_eight(tup)
end

function clip_plane_tet_12(proj::SMatrix{1,4,T,4}, v::NTuple{4,SVector{3,T}}) where {T}
    a = proj_weight_poly(proj, v, 2, 3)
    b = proj_weight_poly(proj, v, 2, 4)
    c = proj_weight_poly(proj, v, 1, 4)
    d = proj_weight_poly(proj, v, 1, 3)
    tup = 0.0 < proj[1] ? (a, b, c, d) : (d, c, b, a)
    return poly_eight(tup)
end

function clip_plane_tet_13(proj::SMatrix{1,4,T,4}, v::NTuple{4,SVector{3,T}}) where {T}
    a = proj_weight_poly(proj, v, 1, 2)
    b = proj_weight_poly(proj, v, 1, 4)
    c = proj_weight_poly(proj, v, 3, 4)
    d = proj_weight_poly(proj, v, 3, 2)
    tup = 0.0 < proj[1] ? (a, b, c, d) : (d, c, b, a)
    return poly_eight(tup)
end

function clip_plane_tet_14(proj::SMatrix{1,4,T,4}, v::NTuple{4,SVector{3,T}}) where {T}
    a = proj_weight_poly(proj, v, 1, 3)
    b = proj_weight_poly(proj, v, 1, 2)
    c = proj_weight_poly(proj, v, 4, 2)
    d = proj_weight_poly(proj, v, 4, 3)
    tup = 0.0 < proj[1] ? (a, b, c, d) : (d, c, b, a)
    return poly_eight(tup)
end
