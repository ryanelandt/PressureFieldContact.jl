
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
