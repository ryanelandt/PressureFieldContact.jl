
function triangle_vertices(i_tri::Int64, m::MeshCache)
    ind_vert = m.tri.ind[i_tri]
    cart_vert = m.point[ind_vert]
    return cart_vert
end

function tetrahedron_vertices_ϵ(i_tet::Int64, m::MeshCache)
    ind_vert = m.tet.tet.ind[i_tet]
    ϵ = m.tet.ϵ[ind_vert]
    cart_vert = m.point[ind_vert]
    return cart_vert, ϵ
end

function calc_ζ_transforms(frame_ζ::CartesianFrame3D, frame_r ::CartesianFrame3D, p_tet, x_r_w, x_w_r)
    x_r_ζ = MatrixTransform(frame_ζ, frame_r, asMatOnePad(p_tet))
    x_w_ζ = x_w_r * x_r_ζ
    x_ζ_w = inv(x_r_ζ) * x_r_w  # NOTE: inv(A_r¹_ζ) is **always** Float64
    return x_w_ζ, x_ζ_w
end

function find_plane_tet(E::Float64, ϵ::SVector{4,Float64}, X_r_w)
    return (E * ϵ) * X_r_w
end
